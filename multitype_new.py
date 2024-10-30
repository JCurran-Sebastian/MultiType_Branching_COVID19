import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from tqdm import tqdm
import numpy.linalg as nla
import seaborn as sns
import scipy as sp
import pandas as pd
import sympy as sym
from odeintw import odeintw


# Helpful function for inverting a matrix - more efficient and avoids numerical error
def invert_mat_safe(A): 
    P, L, U = sp.linalg.lu(A.copy())
    I = np.identity(len(A.copy()))
    Psys = sp.linalg.solve(P, I)
    Lsys = sp.linalg.solve_triangular(L, Psys, lower = True)
    Usys = sp.linalg.solve_triangular(U, Lsys)
    return Usys




def Omega_sym(t):
    a, b, c, d, s, g, beta = sym.symbols('a b c d s g beta')
    Omat = sym.Matrix([[-s, 0, 0, 0, s, 0, 0, 0], 
                     [0, -s, 0, 0, 0, s, 0, 0], 
                     [0, 0, -s, 0, 0, 0, s, 0], 
                     [0, 0, 0, -s, 0, 0, 0, s],
                      [a*beta, b*beta, c*beta, d*beta, -g, 0, 0, 0],
                      [a*beta, b*beta, c*beta, d*beta, 0, -g, 0, 0],
                      [a*beta, b*beta, c*beta, d*beta, 0, 0, -g, 0],
                      [a*beta, b*beta, c*beta, d*beta, 0, 0, 0, -g]]
                     )
    return Omat


def eigvls_sym(params, offspring):
    evals = offspring(0).diagonalize()[1]
    return np.diag(np.array(evals.subs(params)).astype(np.float64))

def eigvects_sym(params, offspring, norm = True):
    evects = offspring(0).diagonalize()[0]
    evects = np.array(evects.subs(params)).astype(np.float64)
    if norm:
        for j in range(0, np.shape(evects)[1]):
            evects[:, j] /= np.sqrt(np.sum(evects[:, j]**2))
    return evects

def reorder_evecs(eigenvalues, eigenvectors, ordering):
    eigenvalues_new = np.asarray([eigenvalues[j] for j in ordering])
    eigenvectors_new = np.asarray([eigenvectors[:,j] for j in ordering]).T
    return eigenvalues_new, eigenvectors_new


def set_odes(u, t, offspring, omega):
    P = offspring
    numeqs = len(omega)
    deriv = np.zeros(numeqs)
    for i in range(0, numeqs):
        deriv[i] = -omega[i]*u[i] + omega[i]*P(u, t)[i]
    return deriv

def set_odes_im(u, t, offspring, omega, eta):
    P = offspring
    numeqs = len(omega)
    deriv = np.zeros(numeqs+1)
    sumqs = 0
    for i in range(0, numeqs):
        deriv[i] = -omega[i]*u[i] + omega[i]*P(u, t)[i]
        sumqs += eta[i]*u[i] 
    deriv[-1] = -np.sum(eta)*u[-1] + u[-1]*sumqs
    return deriv



def set_mean_odes(u, t, Omega):
    deriv = Omega@u
    return deriv
def set_mean_odes_im(u, t, Omega, eta):
    deriv = Omega@u + eta
    return deriv
    

def variance(t, y0, omega, Omega, params, return_vec = True, im = [None]):
    prop_vec, const_vec, sigma, beta_baseline, gamma = params
    if im[0] == None:
        eta = np.zeros(len(omega))
    else: 
        eta = im
    ntypes = len(omega)
    time = t.copy()
    scale = prop_vec*const_vec

    nexposed = int(ntypes/2)
    
    
    eigvls, orth = nla.eig(Omega) 
    if np.iscomplex(orth).any():
        sym_params = [scale[0], scale[1], scale[2], scale[3], sigma, gamma, beta_baseline]
        sym_labels = ['a', 'b', 'c', 'd', 's', 'g', 'beta']
        sym_dict = dict([(l, p) for (l, p) in zip(*(sym_labels, sym_params))])
        diag_sym = Omega_sym(0).diagonalize()
        sym_eigvectors = diag_sym[0]
        sym_eigvalues = diag_sym[1]
        eigvls = np.diag(np.array(sym_eigvalues.subs(sym_dict)))
        orth = np.array(sym_eigvectors.subs(sym_dict))
    assert(not np.iscomplex(orth).any()), 'Jacobian has complex eigendecomposition'
    eigvls = eigvls.astype('float32')
    orth = orth.astype('float32')

    ordering = (np.argsort(eigvls)).tolist()
    eigvls, orth = reorder_evecs(eigvls, orth, ordering)
    diagmat = np.diag(eigvls)
    growth_rate = np.max(eigvls)
    assert(growth_rate > 0), 'Growth Rate must be positive (i.e. Branching Process must be supercritical)'
    
    orth_inv = invert_mat_safe(orth.copy())
    orth_c = orth.copy().T
    orthc_inv = invert_mat_safe(orth_c.copy())
    
    ### Build variance matrix via Kronecker Products
    H = np.kron(orth, np.kron((orthc_inv), orthc_inv))
    Hinv = np.kron(orth_inv, np.kron(orth_c, orth_c)) 
    Amat = np.kron(orthc_inv, orthc_inv)
    Amat_inv = np.kron(orth_c, orth_c)
    
    vec_w = np.zeros(ntypes**3)
    
    # Build Hessian matrix
    Hessian_mat = np.zeros((ntypes, ntypes, ntypes))
    

    for level in range(0, nexposed):
        for i in range(0, nexposed):

            Hessian_mat[level + nexposed, level+nexposed, i] = (beta_baseline * const_vec[i]*prop_vec[i])/(omega[level+nexposed])
            Hessian_mat[level + nexposed, i, level+nexposed] = (beta_baseline * const_vec[i]*prop_vec[i])/(omega[level+nexposed])


    
    
    P_mat = np.zeros((ntypes, ntypes))
    

    for nex in range(0, nexposed):
        P_mat[nex, nex+nexposed] = omega[nex]/omega[nex]
        P_mat[nex+nexposed, :nexposed] = const_vec*prop_vec*beta_baseline/omega[nex+nexposed]
        P_mat[nex+nexposed, nex+nexposed] = np.sum(const_vec * beta_baseline * prop_vec)/omega[nex+nexposed]

    P_mat = P_mat.T
    
    Gmat = np.zeros((ntypes, ntypes, ntypes))
    C = np.zeros(ntypes**3)
    
    for l in range(0, ntypes):

        Gmat[l, :, :] = Hessian_mat[l, :, :] + np.diag(P_mat[:, l]) - np.outer(P_mat[:, l], P_mat[:, l])
        unitvec = np.zeros(ntypes)
        unitvec[l] = 1
        Gmat[l, :, :] += np.outer(unitvec, unitvec) + np.outer(P_mat[:, l], P_mat[:, l]) - np.outer(unitvec, P_mat[:, l]) - np.outer(P_mat[:, l], unitvec)
        Gmat[l, :, :] *= omega[l]
        C[l*ntypes*ntypes:(l+1)*ntypes*ntypes] = Gmat[l, :, :].flatten('F') # Stack columns for 'vec' operator
    
    HinvC = Hinv @ C
    
    # Build diagonal matrices for Kronecker Product
    
    kp = np.kron(diagmat, np.identity(ntypes)) + np.kron(np.identity(ntypes), diagmat) # Kronecker product
    kp_inv = nla.inv(kp)
    
    
    
    diagmat_inv = nla.inv(diagmat)
    idn = np.identity(ntypes)
    dkron = np.kron(-diagmat, np.identity(ntypes**2)) + np.kron(idn, kp)
    dkron_inv = nla.inv(dkron)
    
    var_vec = np.zeros((len(time)))
    var_mat = np.zeros((len(time), ntypes, ntypes))
    
    T_idx = 0
    for T in tqdm(t):
        vec_w = np.zeros(ntypes**3)
        vec_v = np.zeros(ntypes**3)
        
        integral_1 = np.kron(idn, (sp.linalg.expm(T*kp) - sp.linalg.expm(0*kp))@kp_inv)@dkron_inv
        integral_2 = np.kron((sp.linalg.expm(T*diagmat) - idn)@diagmat_inv, np.kron(idn, idn))@dkron_inv
        Dmat_im = integral_1 - integral_2
        Dmat = (sp.linalg.expm(kp*T) - sp.linalg.expm(kp*0))@kp_inv
        vecvar = H @ Dmat_im  @ HinvC
        sum_vecvar_i = np.zeros((ntypes, ntypes))

        for i in range(0, ntypes):
            
            unitvec = np.zeros(ntypes)
            unitvec[i] = 1.
            vecvar_i =  np.kron(unitvec, np.kron(np.identity(ntypes), np.identity(ntypes))) @ vecvar
            sum_vecvar_i += np.real_if_close(y0[i] * np.reshape(vecvar_i, (ntypes, ntypes)).transpose())
            unit_outer_prod = np.outer(unitvec, unitvec).flatten() # Stack columns for 'vec' operator
            vec_wi =  eta[i] * (Amat @ Dmat @ Amat_inv @ unit_outer_prod)
            
            vec_w[i*ntypes*ntypes:(i+1)*ntypes*ntypes] = vec_wi
            vec_v[i*ntypes*ntypes:(i+1)*ntypes*ntypes] = eta[i] * vecvar_i
            
        v = np.reshape(vec_v, (ntypes, ntypes, ntypes)).transpose(0, 2, 1)
        w = np.reshape(vec_w, (ntypes, ntypes, ntypes)).transpose(0, 2, 1)

        v += w
        v += sum_vecvar_i
        
        var_vec[T_idx] = np.sum(v)
        var_mat[T_idx, :, :] = np.sum(v, axis = 0)
        T_idx +=1
    if return_vec:
        return var_vec
    else:
        return var_mat


def Tstar(t, p_zero, var_coeff, thresh1=1e-2, thresh2=1e-2):
    dt = t[1] - t[0]
    diff_pzero = np.abs(np.gradient(p_zero, dt))
    try:
        T1_idx = int(np.min(np.where(diff_pzero<=thresh1)[0]))
    except:
        raise ValueError("Try higher value of thresh1")
    
    diff = diff = np.abs(np.gradient(var_coeff, dt))
    try:
        T2_idx = int(np.min(np.where(diff<=thresh2)[0])) 
    except:
        raise ValueError("Try higher value of thresh2")
    Tstar_idx = int(np.max((T1_idx, T2_idx)))
    Tstar = t[Tstar_idx]
    return [Tstar_idx, Tstar]
