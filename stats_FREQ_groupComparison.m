%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_FREQ.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Statistical analysis of frequency results (between-groups study design) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% = Settings =

% PLEASE SPECIFY the folder for this statistical analysis
stats_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\25vs12\';


% location of freq results for each group:
migraineurs_folder = [stats_folder 'GA_migraineurs\'];
controls_folder = [stats_folder 'GA_controls\'];

% load data
load([migraineurs_folder 'GA_avg.mat']);
mig_avg = GA_freq;
load([controls_folder 'GA_avg.mat']);
ctrl_avg = GA_freq;

load([migraineurs_folder 'GA_individuals.mat']);
mig_indi = GA_freq_indi; 
load([controls_folder 'GA_individuals.mat']);
ctrl_indi = GA_freq_indi;


%% plot overall power (mig vs ctrl)
x_limits = [2 30];

figure; hold on;
plot(mig_avg.freq, mean(mig_avg.powspctrm));
plot(ctrl_avg.freq, mean(ctrl_avg.powspctrm));

xlim(x_limits);
xlabel('Frequency (Hz)');
ylabel('Absolute power (uV^2)');
legend({'Migraineurs', 'Controls'});
hold off;

export_fig(gcf, [stats_folder 'overall_power_mig-vs-ctrl.png']);


%% run t-test (mig vs ctrl) on a particular freq range
freq_range = [9:12];

mig = squeeze(mean(mig_indi.powspctrm, 2)); % avg over all channels
ctrl = squeeze(mean(ctrl_indi.powspctrm, 2)); % avg over all channels

a = mean(mig(:,freq_range), 2); % avg over the selected freq range
b = mean(ctrl(:,freq_range), 2); % avg over the selected freq range

% this function conducts a two-sample t-test
[h,p,ci,stats] = ttest2(a, b, 'Vartype','unequal'); % or should it be equal variance? (p-values are similar)
stats.tstat  % t-value
p            % p-value


%% individual channel analysis: 12vs12 t-test at each freq for each channel

N_chan = size(mig_indi.powspctrm, 2); % number of channels
N_freq = size(mig_indi.powspctrm, 3); % number of freqs

% initialise "chan x freq" matrix to store t-values
t_values = zeros(N_chan, N_freq);

for i = 1:N_chan % loop through each channel
    for j = 1:N_freq % loop through each freq
        a = mig_indi.powspctrm(:,i,j); % extract power for all migraineurs
        b = ctrl_indi.powspctrm(:,i,j); % extract power for all controls
        
        [h,p,ci,stats] = ttest2(a, b, 'Vartype','unequal'); % or should it be equal variance?
        t_values(i,j) = stats.tstat; % store the t-value into the "chan x freq" matrix
    end
end

figure; title('t-values (migraineurs > controls)');
imagesc(t_values)
colorbar
ylabel('EEG channel');
xlabel('Frequency (Hz)');

export_fig(gcf, [stats_folder 'indi-chan-analysis.png']);


%% Figure 2b
freq_range = [9:12];
freq_band = 'alpha';

% for each channel, compute average power over the selected freq range
N_chan = size(mig_indi.powspctrm, 2); % number of channels

% initialise an array to store the t-value for each channel
t_values = zeros(N_chan, 1);

for i = 1:N_chan % loop through each channel
    a = mig_indi.powspctrm(:,i,freq_range); % extract power for all migraineurs (only for the selected freq range)
    a = mean(a,3); % take the mean over that freq range
    b = ctrl_indi.powspctrm(:,i,freq_range); % do the same for controls
    b = mean(b,3); 
    
    [h,p,ci,stats] = ttest2(a, b, 'Vartype','unequal'); % or should it be equal variance?
    t_values(i) = stats.tstat; % store the t-value into the array
end

% create a dummy var for plotting
freq = GA_freq;
freq.powspctrm = t_values;
freq.freq = 0; % we no longer have a frequency dimension, just fill with a dummy value

% plot topography based on the t-values
load('lay_NeuroPrax32.mat');
plot_TFR_topo(freq, lay, freq_band, [], [stats_folder 'indi-chan-analysis_tvalues_'])


%% plot the topography difference btwn two groups (like in Flavia's paper)
% e.g. subtract the topo for alpha (migraineurs minus controls)

% NOTE: this plot is based on absolute power diff in the GA, 
% whereas Flavia plotted the t-values

% calculate the difference GA
diff_GA = mig_avg;
diff_GA.powspctrm = mig_avg.powspctrm - ctrl_avg.powspctrm;

% plot topography based on the difference GA
load('lay_NeuroPrax32.mat');
plot_TFR_topo(diff_GA, lay, 'theta', [4 8], [stats_folder 'topo_mig-minus-ctrl_'])
plot_TFR_topo(diff_GA, lay, 'alpha', [9 12], [stats_folder 'topo_mig-minus-ctrl_'])
plot_TFR_topo(diff_GA, lay, 'beta', [13 25], [stats_folder 'topo_mig-minus-ctrl_'])
