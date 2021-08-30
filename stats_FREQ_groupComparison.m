% location of results folder for each group:
migraineurs_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\12vs12\GA_migraineurs\';
controls_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\12vs12\GA_controls\';


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


%% run t-test (mig vs ctrl) on a particular freq range
freq_range = [8:10];

load([migraineurs_folder 'GA_individuals.mat']);
mig = squeeze(mean(GA_freq_indi.powspctrm, 2)); % avg over all channels
load([controls_folder 'GA_individuals.mat']);
ctrl = squeeze(mean(GA_freq_indi.powspctrm, 2)); % avg over all channels

a = mean(mig(:,freq_range), 2); % avg over the selected freq range
b = mean(ctrl(:,freq_range), 2); % avg over the selected freq range
[h,p,ci,stats] = ttest2(a, b, 'Vartype','unequal'); % or should it be equal variance? (p-values are similar)
p
