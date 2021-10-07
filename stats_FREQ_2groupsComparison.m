%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_FREQ_2groupsComparison.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Statistical analysis of frequency results:
% Comparison between 2 groups (e.g. patients vs controls)
%
% TODO: can merge this script into stats_FREQ_nGroups.m
% (note - 2 groups is a special case, need to use t-test rather than F-test)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% = Settings =

% Run t-tests using equal variance (i.e. pooled) or unequal variance? (p-values obtained were similar)
% 'unequal' is the most conservative approach
% https://www.investopedia.com/terms/t/t-test.asp
% https://statmagic.info/Content/Help-Content/two-sample-mean.html
varType = 'unequal'; %'equal';

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

load('lay_NeuroPrax32.mat');


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
[h,p,ci,stats] = ttest2(a, b, 'Vartype',varType); % or should it be equal variance? (p-values were similar)
stats.tstat  % t-value
p            % p-value


%% individual channel analysis: 25vs12 t-test at each freq for each channel
% (Figure 2a in Flavia paper)

N_chan = size(mig_indi.powspctrm, 2); % number of channels
N_freq = size(mig_indi.powspctrm, 3); % number of freqs

% initialise "chan x freq" matrix to store t-values
t_values = zeros(N_chan, N_freq);

for i = 1:N_chan % loop through each channel
    for j = 1:N_freq % loop through each freq
        a = mig_indi.powspctrm(:,i,j); % extract power for all migraineurs
        b = ctrl_indi.powspctrm(:,i,j); % extract power for all controls
        
        [h,p,ci,stats] = ttest2(a, b, 'Vartype',varType);
        t_values(i,j) = stats.tstat; % store the t-value into the "chan x freq" matrix
    end
end

figure; title('t-values (migraineurs > controls)');
imagesc(t_values)
colorbar
ylabel('EEG channel');
xlabel('Frequency (Hz)');

export_fig(gcf, [stats_folder 'indi-chan-analysis\indi-chan-analysis.png']);


%% Figure 2b in Flavia paper
freq_range = [9:12];
freq_band = 'alpha';

% for each channel, compute average power over the selected freq range
N_chan = size(mig_indi.powspctrm, 2); % number of channels

% initialise an array to store the t-value for each channel
t_values = zeros(N_chan, 1);
p_values = zeros(N_chan, 1);

for i = 1:N_chan % loop through each channel
    a = mig_indi.powspctrm(:,i,freq_range); % extract power for all migraineurs (only for the selected freq range)
    a = mean(a,3); % take the mean over that freq range
    b = ctrl_indi.powspctrm(:,i,freq_range); % do the same for controls
    b = mean(b,3); 
    
    [h,p,ci,stats] = ttest2(a, b, 'Vartype',varType); % or should it be equal variance?
    t_values(i) = stats.tstat; % store the t-value into the array
    p_values(i) = p;
end

% create a dummy var for plotting
freq = GA_freq;
freq.powspctrm = t_values;
freq.freq = 0; % we no longer have a frequency dimension, just fill with a dummy value

% plot topography based on the t-values
plot_TFR_topo(freq, lay, freq_band, [], [stats_folder 'indi-chan-analysis\tvalues_']);


%% plot the topography difference btwn two groups
% e.g. subtract the topo for alpha (migraineurs minus controls)

% NOTE: this plot is based on absolute power diff in the GA, 
% whereas Flavia plotted the t-values (see above - "Figure 2b")

% calculate the difference GA
diff_GA = mig_avg;
diff_GA.powspctrm = mig_avg.powspctrm - ctrl_avg.powspctrm;

% plot topography based on the difference GA
plot_TFR_topo(diff_GA, lay, 'theta', [4 8], [stats_folder 'topo_mig-minus-ctrl_'])
plot_TFR_topo(diff_GA, lay, 'alpha', [9 12], [stats_folder 'topo_mig-minus-ctrl_'])
plot_TFR_topo(diff_GA, lay, 'beta', [13 25], [stats_folder 'topo_mig-minus-ctrl_'])



%% Cluster-based statistical analysis
% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/
% https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_stats/#2-compute-between-participants-contrasts

% load the data
load([migraineurs_folder 'allSubjects_freq.mat']);
mig = allSubjects_freq;
load([controls_folder 'allSubjects_freq.mat']);
ctrl = allSubjects_freq;

load('neighbours_NeuroPrax32.mat'); % obtained using 'trigangulation' method in ft_prepare_neighbour


% Opt 1: find spatio-freq cluster:
%foi = [1 30];
%freq_band = 'spatio-freq';
% Opt 2: find spatial cluster for a particular freq band:
foi = 9; %[9 12];
freq_band = 'alpha';

cfg = [];
cfg.channel = 'all';
cfg.frequency = foi;
cfg.avgoverfreq = 'yes'; % enable for Opt 2

cfg.method = 'montecarlo';
cfg.correctm = 'cluster'; %'no'; % it is common in MEG studies to run uncorrected at cfg.alpha = 0.001
cfg.clusteralpha = 0.05; % threshold for selecting candidate samples to form clusters
cfg.clusterstatistic = 'maxsum';
%cfg.clusterthreshold = 'nonparametric_common';
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment
cfg.minnbchan = 2; % minimum number of neighbourhood channels required to be significant 
                   % in order to form a cluster 
                   % (default: 0, ie. each single channel can be considered a cluster)

cfg.alpha = 0.1; %0.001  % threshold for cluster-level statistics (any cluster with a p-value lower than this will be reported as sig - an entry of '1' in .mask field)
cfg.numrandomization = 2000; % Rule of thumb: use 500, and double this number if it turns out 
    % that the p-value differs from the chosen alpha (e.g. 0.05) by less than 0.02
    
%{
% within-subject design
numSubjects = length(data.(eventnames_real{1})); % check how many subjects we are including
within_design_1x2 = zeros(2,2*numSubjects);
within_design_1x2(1,:) = repmat(1:numSubjects,1,2);
within_design_1x2(2,1:numSubjects) = 1;
within_design_1x2(2,numSubjects+1:2*numSubjects) = 2;

within_design_1x3 = zeros(2, 3*numSubjects);
within_design_1x3(1, :) = repmat(1:numSubjects, 1, 3);
within_design_1x3(2, 1:numSubjects) = 1;
within_design_1x3(2, numSubjects+1:2*numSubjects) = 2;
within_design_1x3(2, 2*numSubjects+1:3*numSubjects) = 3;

cfg.uvar  = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar  = 2; % row of design matrix that contains independent variable (i.e. the conditions)
%}

% design matrix
num_mig = length(mig);
num_ctrl = length(ctrl);

design = zeros(1, num_mig + num_ctrl);
design(1, 1:num_mig) = 1;
design(1, (num_mig+1):(num_mig + num_ctrl)) = 2;

cfg.design = design;
cfg.ivar   = 1;

% Make sure we are using 2-tailed t-tests:
cfg.statistic = 'indepsamplesT'; % independent samples t-test (patients vs controls)
cfg.tail = 0;
cfg.clustertail = 0; % 2 tailed test
cfg.correcttail = 'prob'; % correct for 2-tailedness

[stat] = ft_freqstatistics(cfg, mig{:}, ctrl{:});
length(find(stat.mask)) % display how many chans were significant/marginal

%save([stats_folder 'cluster_stat\minnbchan' mat2str(cfg.minnbchan) '_' freq_band '.mat'], 'stat');


%% ft_clusterplot (plots t-values by default)
% this is a wrapper around ft_topoplot, automatically extracts info about the cluster

load([stats_folder 'cluster_stat\minnbchan2_alpha.mat']);

% use a nice-looking colourmap
%ft_hastoolbox('brewermap', 1); % ensure this toolbox is on the path
%cmap = colormap(flipud(brewermap(64, 'RdBu')));

% too much warning, can't see the console output
ft_warning off 'FieldTrip:ft_clusterplot:ft_topoplotTFR:topoplot_common:ft_selectdata:getdimord:warning_dimord_could_not_be_determined:line681'

cfg = [];
%cfg.zlim = [-5 5]; % set scaling (range of t-values) (usually using automatic is ok) 
cfg.highlightcolorpos = [0 0 0]; % black for pos clusters
cfg.highlightcolorneg = [255/255 192/255 203/255]; % pink for neg clusters
cfg.alpha = 0.1; % any clusters with a p-value below this threshold will be plotted
cfg.layout = lay;
cfg.style = 'straight';
%cfg.colormap = cmap;

% turn on the following lines if you are after one particular subplot
cfg.subplotsize = [1 1];
cfg.highlightsizeseries = [9 9 9 9 9]; % make the highlight markers bigger
cfg.colorbar = 'yes'; % shows the scaling

ft_clusterplot(cfg, stat);

%export_fig(gcf, [stats_folder 'indi-chan-analysis\cluster_stat\minnbchan2_' freq_band '.png']);
