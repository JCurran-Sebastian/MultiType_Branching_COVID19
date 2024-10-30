import numpy as np
# import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy.linalg as nla
import scipy as sp
import pandas as pd
import sympy as sym
import multitype_new as mt
import importlib
importlib.reload(mt)
from scipy.integrate import odeint

# Parameters to change 
# Change idx to change value of R from Reff_vec - or can choose a different value of R altogether by changing the value of RV

ntypes = 8
# Number of types-at-birth
nexposed = int(ntypes/2)


# Disease Parameters

gamma = 0.25 #0.4 # Recovery Rate
sigma = 1/3 #0.3 # Progression rate from exposed to infectious

beta_baseline = np.squeeze(pd.read_csv('VOC_beta.csv',header=None).to_numpy())

print('Running for beta = ' + str(beta_baseline))

epsilon = 1e-3
RV = beta_baseline/gamma

# Relative susceptibility for different types 
#sus_ur = 1-0.4 # Susceptibility of unvaccinated, previously infected against new strain
#sus_vu = 1-0.6 # Susceptibility of vaccinated, not previously infected against new strain
#sus_vr = np.min((sus_ur, sus_vu)) # Susceptibility of vaccinated, previously infected against new strain


Tmax = 1000 
num = Tmax 
time = np.linspace(0, Tmax, num=num)
dtime = time[1]-time[0]


# Proportion of population in each type
# p_vac = 0.7
# p_no_vac = 1- p_vac 

# p_rec = 0.2 # % recovered from previous strains
# p_sus = 1-p_rec # % never infected with any coronavirus

# Reduced susceptibility for each type based on infection/vaccine immunity
const_vec = np.squeeze(pd.read_csv('const_vec_in.csv',header=None).to_numpy()) #np.array((1., sus_ur, sus_vu, sus_vr))
prop_vec = np.squeeze(pd.read_csv('prop_vec_in.csv',header=None).to_numpy()) #np.array((p_sus*p_no_vac, p_rec*p_no_vac, p_sus*p_vac, p_rec*p_vac))
scale = const_vec*prop_vec


# Vector of immigrants at each type - currently only 1 immigrant of type 1 (unvaccinated, no prior infection) every 5 days

im = 0
im_vec = np.zeros(ntypes)
y0 = np.zeros(ntypes)
y0[0] = 1
ics = np.zeros(ntypes)
def eta(t):
    return im_vec

# No reduction in transmission conferred by vaccination
trans_vec = np.ones(nexposed)

def lifetime_vec(beta, gamma, sigma, ntypes = ntypes, const_vec = const_vec, prop_vec = prop_vec, trans_vec = trans_vec):
    nexposed = int(ntypes / 2)
    omega_vec = sigma * np.ones(ntypes) 

        
    for nex in range(0, nexposed):
        omega_vec[nex + nexposed] = np.sum(const_vec*prop_vec*beta)*trans_vec[nex] + gamma
    return omega_vec

omega_vec = lifetime_vec(beta_baseline, gamma, sigma)

def P(u, t): # Offspring distribution

        ## Returns a vector of length ntypes with entry i containing  probabilities of generating particles 
        ## of each type from a particle of type i

        pvec = np.zeros_like(u)
        nexposed = int(ntypes/2)
        omega = (np.sum(prop_vec * const_vec * beta_baseline) + gamma)
        

        sum_gen_funcs = np.sum((beta_baseline*prop_vec/omega)*const_vec*u[:nexposed])
         
        for nex in range(0, nexposed):

            pvec[nex] = u[nex+nexposed] 
            pvec[nex+nexposed] =  u[nex+nexposed] * sum_gen_funcs  + (gamma/omega)
        
        return pvec  

def Jacobian_mat(beta, gamma, sigma, ntypes = ntypes, const_vec = const_vec, prop_vec = prop_vec, trans_vec = trans_vec, omega_vec = omega_vec):
    Omat = np.zeros((ntypes, ntypes))
    omega_vec = sigma * np.ones(ntypes)       
    nexposed = int(ntypes/2)
    for nex in range(0, nexposed):
        omega_vec[nex + nexposed] = np.sum(const_vec*prop_vec*beta)*trans_vec[nex] + gamma
        Omat[nex, nex+nexposed] = omega_vec[nex]
        Omat[nex+nexposed, :nexposed] = const_vec*prop_vec*beta *trans_vec[nex]
        Omat[nex+nexposed, nex+nexposed] = np.sum(const_vec * beta * prop_vec)*trans_vec[nex]

    Omat -= np.diag(omega_vec)
    return Omat
	
	
rho = (odeint(mt.set_odes, ics, time, args = (P, omega_vec)).T )
q = np.ones_like(time)
for i in range(ntypes):
    q *= rho[i, :] ** y0[i]

np.savetxt('./Outputs_for_matlab/p_extinction/p_extinction_R=' + str(RV) + '.csv', np.vstack((time, q)))

# Jacobian, or mean matrix
Omat = Jacobian_mat(beta_baseline, gamma, sigma)
# Mean of each type
Mmat = odeint(mt.set_mean_odes, y0, time, args = (Omat,)).T
# Total mean
mean = np.sum(Mmat, axis=0)

# print(y0, omega_vec)
# print(pd.DataFrame(Omat))
# print([prop_vec, const_vec, beta_baseline, gamma])
# print(eta(10))


variance = mt.variance(time, y0, omega_vec, Omat, [prop_vec, const_vec, sigma, 
                                                   beta_baseline, gamma], eta)

#Coefficient of variation
sig_over_mean = (np.sqrt(variance)/(mean))
np.savetxt('./Outputs_for_matlab/mean_variance/mean_variance_R=' + str(RV) + '.csv', np.vstack((time, mean, variance)))
	
Tstar_idx, Tstar = mt.Tstar(time, q, sig_over_mean, thresh1 = epsilon, thresh2 = epsilon)
Zstar_min= mean[Tstar_idx]

eigvls, orth = nla.eig(Omat)

if np.iscomplex(orth).any():

# Maybe change this back
    print('Numerical instability: Switching to symbolic compution of eigenvectors')
    sym_params = [scale[0], scale[1], scale[2], scale[3], sigma, gamma, beta_baseline]
    sym_labels = ['a', 'b', 'c', 'd', 's', 'g', 'beta']
    sym_dict = dict([(l, p) for (l, p) in zip(*(sym_labels, sym_params))])
    diag_sym = mt.Omega_sym(0).diagonalize()
    sym_eigvectors = diag_sym[0]
    sym_eigvalues = diag_sym[1]
    eigvls = np.diag(np.array(sym_eigvalues.subs(sym_dict)))
    orth = np.array(sym_eigvectors.subs(sym_dict))
    for col in range(orth.shape[0]):
        orth[:, col] /= np.sum(orth[:, col]**2)
assert(not np.iscomplex(orth).any()), 'Jacobian has complex eigendecomposition'
eigvls = eigvls.astype('float64')
orth = orth.astype('float64')

# Some issue with ordering in sympy vs numpy
ordering = np.flip(np.argsort(eigvls)).tolist()
eigvls, orth = mt.reorder_evecs(eigvls, orth, ordering)


diagmat = np.diag(eigvls)
growth_rate = np.max(eigvls)
assert(growth_rate > 0), 'Growth Rate must be positive (i.e. Branching Process must be supercritical)'

orth_inv = mt.invert_mat_safe(orth.copy())
orth_c = orth.copy().T
orthc_inv = mt.invert_mat_safe(orth_c.copy())


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

        Hessian_mat[level + nexposed, level+nexposed, i] = (beta_baseline * const_vec[i]*prop_vec[i])/(omega_vec[level+nexposed])
        Hessian_mat[level + nexposed, i, level+nexposed] = (beta_baseline * const_vec[i]*prop_vec[i])/(omega_vec[level+nexposed])

P_mat = np.zeros((ntypes, ntypes))


for nex in range(0, nexposed):
    P_mat[nex, nex+nexposed] = omega_vec[nex]/omega_vec[nex]
    P_mat[nex+nexposed, :nexposed] = const_vec*prop_vec*beta_baseline/omega_vec[nex+nexposed]
    P_mat[nex+nexposed, nex+nexposed] = np.sum(const_vec * beta_baseline * prop_vec)/omega_vec[nex+nexposed]

P_mat = P_mat.T

Gmat = np.zeros((ntypes, ntypes, ntypes))
C = np.zeros(ntypes**3)
blorg = np.zeros((ntypes, ntypes))
for l in range(0, ntypes):
    f_vec = (P_mat - np.identity(ntypes))[l, :]
    unitvec = np.zeros(ntypes)
    unitvec[l] = 1

    Gmat[l, :, :] = Hessian_mat[l, :, :] - np.diag(P_mat[:, l]) #- np.outer(P_mat[:, l], P_mat[:, l])
    Gmat[l, :, :] += np.diag(unitvec)        
    Gmat[l, :, :] *= omega_vec[l]


growth_rate = np.max(eigvls)
print(growth_rate)

#####
##### You have changed the minus sign here!!!! Change back if things go wrong.
##### Update: You have changed it back. If in doubt, it should have a minus sign in front
#####

if (orth[:, 0] < 0).all():
    orth = -orth
change_from_ebasis = orth # P

evec = change_from_ebasis[:, 0]

change_to_ebasis = mt.invert_mat_safe(change_from_ebasis)
xi = change_to_ebasis[:, 0]
variance_diffusion_mat = np.zeros((ntypes, ntypes))

for l in range(ntypes): 
    variance_diffusion_mat += Gmat[l, :, :] *evec[l]  
variance_all = (orth[:, 0].T@variance_diffusion_mat@orth[:, 0])

print('Total variance is ' + str(variance_all))

evec_scaling = np.abs(1/np.sum(evec))

# change_from_ebasis = -orth # P
# evec = change_from_ebasis[:, 0]
# if (evec<=0).all():
#     evec = np.abs(evec)
# change_to_ebasis = mt.invert_mat_safe(change_from_ebasis)
# xi = change_to_ebasis[:, 0]
# variance_diffusion_mat = np.zeros((ntypes, ntypes))
# for l in range(ntypes):
#     variance_diffusion_mat += change_from_ebasis@(change_from_ebasis)@(Gmat[l, :, :] * xi[l])
    

# variance_all = nla.norm(variance_diffusion_mat[0, :])
# print(variance_all)
# evec_scaling = 1/np.sum(xi)

xvec = np.linspace(np.finfo(float).eps, 20000, int(1e5))
xvec_idx = np.min(np.where(xvec>=(evec_scaling*Zstar_min))[0])
integral_limit = xvec[xvec_idx]

## Save growth rate, variance and scaling factor
print('Zstar = ' + str(integral_limit))

np.savetxt('./Outputs_for_matlab/FPT_params_beta='+str(beta_baseline)+'.csv', np.array((growth_rate, variance_all, evec_scaling, integral_limit,beta_baseline)))
np.savetxt('./Outputs_for_matlab/dominant_eigenvector_beta='+str(beta_baseline)+'.csv', evec)
np.savetxt('./Outputs_for_matlab/Jacobian_beta='+str(beta_baseline)+'.csv', pd.DataFrame(Omat))