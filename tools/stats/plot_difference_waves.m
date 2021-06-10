% Run this after obtaining the relevant stats output, e.g. SwCost_interaction
% (line 240 in stats_ROI.m)


% Plot 3 difference waves
cfg = [];
cfg.channel   = {'all'}; % there is only one channel (i.e. the virtual sensor for this ROI)
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'no'; % average across subjects
GA_nat = ft_timelockgrandaverage(cfg, timelock_SwCost_Nat{:}); 
GA_art = ft_timelockgrandaverage(cfg, timelock_SwCost_Art{:}); 
GA_bi = ft_timelockgrandaverage(cfg, timelock_SwCost_Bi{:}); 

figure('Name', ['GA in ' ROI_name]); hold on
for j = 1:length(eventnames_real)
    plot(GA_nat.time, GA_nat.avg, 'g', 'lineWidth',2);
    plot(GA_art.time, GA_art.avg, 'b', 'lineWidth',2);
    plot(GA_bi.time, GA_bi.avg, 'r', 'lineWidth',2); 
end
xlim(PLOT_XLIM);
legend({'Nat switch cost', 'Art switch cost', 'Bi switch cost'}, 'Location','northwest');


% Plot t-values
figure; hold on;
x = SwCost_interaction.RIFG.time;
plot(x, SwCost_interaction.RIFG.stat, 'k', 'lineWidth',2); 
plot(x, repmat(2, [1,length(x)]), 'k--'); % make a horizontal line at t=2
xlim([0 0.6])
    