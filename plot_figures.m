clear
close all


UK_pop = 66000000;
time_vec = linspace(0, 1000, 1000);
VOC_introduction_date = 200;
datenum(2021,5,17)+VOC_introduction_date + time_vec;
beta_high = 1.1;
beta_low = 0.675;

%% Read in hybrid realisations from running get_new_model_results.m
% Low values of beta
beta_curves_low = readmatrix(append("Plotting_Outputs/save_curves_beta=", string(beta_low), ".csv"));
times_low = load(append("Plotting_Outputs/save_dates_beta=", string(beta_low), ".mat")).save_dates;
%times_reached_up_low = readmatrix(append('Plotting_Outputs/time_reached_up_beta=', string(beta_low), '.csv'));
%times_reached_down_low = readmatrix(append('Plotting_Outputs/time_reached_down_beta=', string(beta_low), '.csv'));
I_hybrid_low = readmatrix(append("Plotting_Outputs/median_hybrid_curve_beta=", string(beta_low), ".csv"));
I_hybrid_times_low = load(append("Plotting_Outputs/median_hybrid_times_beta=", string(beta_low), ".mat")).median_dates_VOC;
I_VOC_low = readmatrix(append('Plotting_Outputs/I_VOC_beta=', string(beta_low),'.csv'));

% High values of beta
beta_curves_high = readmatrix(append("Plotting_Outputs/save_curves_beta=", string(beta_high), ".csv"));
times_high = load(append("Plotting_Outputs/save_dates_beta=", string(beta_high), ".mat")).save_dates;
%times_reached_up_high = readmatrix(append('Plotting_Outputs/time_reached_up_beta=', string(beta_high), '.csv'));
%times_reached_down_high = readmatrix(append('Plotting_Outputs/time_reached_down_beta=', string(beta_high), '.csv'));
I_hybrid_high = readmatrix(append("Plotting_Outputs/median_hybrid_curve_beta=", string(beta_high), ".csv"));
I_hybrid_times_high = load(append("Plotting_Outputs/median_hybrid_times_beta=", string(beta_high), ".mat")).median_dates_VOC;
I_VOC_high = readmatrix(append('Plotting_Outputs/I_VOC_beta=', string(beta_high),'.csv'));



% Read in deterministic outputs
I_UK_dates = load('Plotting_Outputs/I_VOC_dates.mat').I_VOC_dates;
I_UK = readmatrix('Plotting_Outputs/I_UK.csv');
I_VOC_dates = load('Plotting_Outputs/I_VOC_dates.mat').I_VOC_dates;


% Read in First-Passage simulation outputs and analytic Feller approximation

simulations_low = readmatrix(append('Plotting_Outputs/gillespie_sims_beta=', string(beta_low),'.csv'));
feller_pdf_low = readmatrix(append('Plotting_Outputs/gillespie_pdf_beta=', string(beta_low),'.csv'));
simulations_high = readmatrix(append('Plotting_Outputs/gillespie_sims_beta=', string(beta_high),'.csv'));
feller_pdf_high = readmatrix(append('Plotting_Outputs/gillespie_pdf_beta=', string(beta_high),'.csv'));


max_plot_y = max(I_VOC_high) + 0.2 * 1e6;
max_plot_y_sims = max(feller_pdf_high) + 0.005;
to_plot = 1:1000;
start_date = I_UK_dates(16); %1st June 2021
end_date = I_UK_dates(850); %1st July 2023
end_date_fpt = I_UK_dates(289); %1st March 2022


npoints = 2000;
%nsims = length(first_passage_times);
horz_points_low = linspace(0,max(beta_curves_low(:)),npoints+2);
horz_points_low = horz_points_low(2:end-1);

horz_points_high = linspace(0,max(beta_curves_high(:)),npoints+2);
horz_points_high = horz_points_high(2:end-1);

f = figure(gcf);
%f.WindowState = 'maximized';
hold on;
box on;
subplot(2,2,1);
p1 = histogram(simulations_low, 100, 'Normalization','pdf');
hold on
p2 = plot(time_vec, feller_pdf_low, 'Color', [0.6350 0.0780 0.1840,1], 'LineWidth', 2);
legend([p1,p2],'Simulation','Analytic')

ylim([0, max_plot_y_sims])
xlim([0, 250])
xlabel('Days after VOC introduction')
ylabel('Density')
%xlim([start_date, end_date_fpt])
title(append('First-Passage Time to Z^*, \beta =', string(beta_low)))

subplot(2,2,2);
 %rand(1,1000)<0.3;
p2 = plot(times_low(VOC_introduction_date:end, :),beta_curves_low(VOC_introduction_date:end, :),'Color', [0.6350 0.0780 0.1840, 0.01], 'LineWidth', 2);
hold on
p2 = plot(times_high,beta_curves_high*NaN,'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1);
%p3 = plot(median(times_low,2),beta_curves_low(:,1),'Color',"#77AC30",'LineWidth',2);
%p3 = plot(median(times_reached_up_low,2,'omitnan'),horz_points_low*UK_pop,'Color',"#77AC30",'LineWidth',1);

p4 = plot(I_VOC_dates(VOC_introduction_date:end), I_VOC_low(VOC_introduction_date:end),'Color',"#0072BD",'LineWidth',2);
p3=plot(I_hybrid_times_low(VOC_introduction_date:end),I_hybrid_low(VOC_introduction_date:end)*UK_pop, 'LineStyle', '--', 'Color',"#77AC30",'LineWidth',2);
p1 = plot(I_UK_dates(1:300), I_UK(1:300),'Color',"#7E2F8E",'LineWidth',2);
p5 = xline(I_UK_dates(VOC_introduction_date), '--', 'VOC Introduction Date');
legend([p1,p2(end),p3, p4],'Deterministic (Resident)','Hybrid realisations','Median hybrid', 'Deterministic (VOC)')
%legend([p1,p2(1), p3],'Deterministic','Hybrid realisations','Median hybrid')
ylabel('Infections')
xlabel('Time')
ylim([0, max_plot_y])
xlim([start_date, end_date])
title(append('VOC waves, \beta = ', string(beta_low)))



subplot(2,2,3);
p1 = histogram(simulations_high, 'Normalization','pdf');
hold on
p2=plot(time_vec, feller_pdf_high, 'Color', [0.6350 0.0780 0.1840,1], 'LineWidth', 2);
legend([p1,p2],'Simulation','Analytic')
xlim([0, 250])
ylim([0, max_plot_y_sims])
xlabel('Days after VOC introduction')
ylabel('Density')
title(append('First-Passage Time to Z^*, \beta =', string(beta_high)))


subplot(2,2,4);
to_plot = 1:1000; %rand(1,1000)<0.3;
p2 = plot(times_high(VOC_introduction_date:end, :),beta_curves_high(VOC_introduction_date:end, :),'Color', [0.6350 0.0780 0.1840, 0.01], 'LineWidth', 5);
hold on
p2 = plot(times_high,beta_curves_high*NaN,'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1);

%p3 = plot(median(times_reached_up_high,2,'omitnan'),horz_points_high*UK_pop,'Color',"#77AC30",'LineWidth',1);
%plot(median(times_reached_down_high,2,'omitnan'),horz_points_high*UK_pop,'Color',"#77AC30",'LineWidth',1);
p4 = plot(I_VOC_dates(VOC_introduction_date:end), I_VOC_high(VOC_introduction_date:end),'Color',"#0072BD",'LineWidth',2);
p3 = plot(I_hybrid_times_high(VOC_introduction_date:end),I_hybrid_high(VOC_introduction_date:end)*UK_pop, 'LineStyle', '--', 'Color',"#77AC30",'LineWidth',2);
p1 = plot(I_UK_dates(1:300), I_UK(1:300),'Color',"#7E2F8E",'LineWidth',2);
p5 = xline(I_UK_dates(VOC_introduction_date), '--', 'VOC Introduction Date');
legend([p1,p2(end),p3, p4],'Deterministic (Resident)','Hybrid realisations','Median hybrid', 'Deterministic (VOC)')
%legend([p1,p2(1), p3],'Deterministic','Hybrid realisations','Median hybrid')
ylabel('Infections')
xlabel('Time')
ylim([0, max_plot_y])
xlim([start_date, end_date])
title(append('VOC waves, \beta = ', string(beta_high)))


pos = get(gcf, 'Position');
set(gcf,'Position',pos+[-500 -500 500 500]);
%exportgraphics(gcf,'figure_2.png')









