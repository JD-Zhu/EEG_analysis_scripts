% SPECIFY the folder for this statistical analysis
stats_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\12vs12\';

% location of freq results for each group:
migraineurs_folder = [stats_folder 'GA_migraineurs\'];
controls_folder = [stats_folder 'GA_controls\'];


%% plot overall power (mig vs ctrl)
x_limits = [2 30];

figure; hold on;

load([migraineurs_folder 'GA_avg.mat']);
plot(GA_freq.freq, mean(GA_freq.powspctrm));
load([controls_folder 'GA_avg.mat']);
plot(GA_freq.freq, mean(GA_freq.powspctrm));

xlim(x_limits);
xlabel('Frequency (Hz)');
ylabel('Absolute power (uV^2)');
legend({'Migraineurs', 'Controls'});
hold off;

export_fig(gcf, [stats_folder 'overall_power_mig-vs-ctrl.png']);


%% run t-test (mig vs ctrl) on a particular freq range
freq_range = [8:10];

load([migraineurs_folder 'GA_individuals.mat']);
mig = squeeze(mean(GA_freq_indi.powspctrm, 2)); % avg over all channels
load([controls_folder 'GA_individuals.mat']);
ctrl = squeeze(mean(GA_freq_indi.powspctrm, 2)); % avg over all channels

a = mean(mig(:,freq_range), 2); % avg over the selected freq range
b = mean(ctrl(:,freq_range), 2); % avg over the selected freq range

% this function conducts a two-samples t-test
[h,p,ci,stats] = ttest2(a, b, 'Vartype','unequal'); % or should it be equal variance? (p-values are similar)
p


%% plot the topography difference btwn two groups (like in Flavia's paper)
% e.g. subtract the topo for alpha (migraineurs minus controls)
load([migraineurs_folder 'GA_avg.mat']);
mig_avg = GA_freq;
load([controls_folder 'GA_avg.mat']);
ctrl_avg = GA_freq;

% calculate the difference GA
diff = mig_avg;
diff.powspctrm = mig_avg.powspctrm - ctrl_avg.powspctrm;

% plot topography based on the difference GA
load('lay_NeuroPrax32.mat');
plot_TFR_topo(diff, lay, 'theta', [4 8], [stats_folder 'topo_mig-minus-ctrl_'])
plot_TFR_topo(diff, lay, 'alpha', [9 12], [stats_folder 'topo_mig-minus-ctrl_'])
plot_TFR_topo(diff, lay, 'beta', [13 25], [stats_folder 'topo_mig-minus-ctrl_'])
