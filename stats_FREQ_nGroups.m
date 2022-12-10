%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_FREQ_nGroups.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Statistical analysis of frequency results:
% Comparison across 3 groups (e.g. migraine phases / frequency of attacks)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global LAYOUT_FILE; global NEIGHBOURS_FILE; global PLOT_XLIM;
global FREQ_BANDS;
common();

load(LAYOUT_FILE);


% = Settings =

% Use logged power or absolute power?
logged = true;

% PLEASE SPECIFY the folder for this statistical analysis
stats_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\migraine_phases\';
%stats_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\migraine_frequency\';

% PLEASE SPECIFY the subject groups for comparison
groups = {'GA_prodrome', 'GA_postdrome+1', 'GA_interictal'};
%groups = {'GA_lessThan1day', 'GA_1-2days', 'GA_moreThan3days'};
% note: you need to have these GA folders ready inside the stats_folder

% Also set up the figure legends (MAKE SURE the order here is same as in the "groups" variable above)
figure_legends = {'Prodrome', 'Postdrome', 'Interictal'};
%figure_legends = {'< 1 day', '1-2 days', '> 3 days'};


% do not change the stuff below
if logged
    logged_suffix = '_logged';
else
    logged_suffix = '';
end

% locations to save results
stats_folder_indi_chan = [stats_folder '\indi-chan-analysis\'];
mkdir(stats_folder_indi_chan);
stats_folder_cluster = [stats_folder_indi_chan '\cluster_stat' logged_suffix '\'];
mkdir(stats_folder_cluster);

    
%% plot overall power (one line == one group)

% plot the GA for each group
figure; hold on;
for g = 1:length(groups)
    load([stats_folder groups{g} '\GA_avg.mat']); % GA for this group
    plot(GA_freq.freq, mean(GA_freq.powspctrm));
end

xlim(PLOT_XLIM);
xlabel('Frequency (Hz)');
ylabel('Absolute power (uV^2)');
legend(figure_legends);
hold off;

export_fig(gcf, [stats_folder 'overall_power_for_each_group.png']);


% plot log-transformed version
figure; hold on;
for g = 1:length(groups)
    load([stats_folder groups{g} '\GA_avg.mat']); % GA for this group
    GA_freq.powspctrm = log(GA_freq.powspctrm);
    plot(GA_freq.freq, mean(GA_freq.powspctrm));
end

xlim(PLOT_XLIM);
xlabel('Frequency (Hz)');
ylabel('Power (log[uV^2]');
legend(figure_legends);
hold off;

export_fig(gcf, [stats_folder 'overall_power_(logged)_for_each_group.png']);


%% Analysis of overall power
% run ANOVA (comparison across all groups) on a particular freq range

% loop thru each freq band
for band = 1:length(FREQ_BANDS)
    freq_band = cell2mat(FREQ_BANDS{band, 1}); % first field is the freq band name
    freq_range = FREQ_BANDS{band, 2}; % second field is the freq range in Hz
    
    % find the start index & end index for the selected freq range
    load([stats_folder groups{1} '\GA_individuals.mat']); % load one group as an example
    freq_min_idx = find(GA_freq_indi.freq == freq_range(1));
    freq_max_idx = find(GA_freq_indi.freq == freq_range(end));

    
    % collate data for anova (each subject should have a single value)
    data_for_anova = [];
    grouping_var = {};
    for g = 1:length(groups)
        load([stats_folder groups{g} '\GA_individuals.mat']); % GA for this group
        if logged % apply log transformation if needed
            GA_freq_indi.powspctrm = log(GA_freq_indi.powspctrm);
        end
        data = squeeze(mean(GA_freq_indi.powspctrm, 2)); % avg over all channels
        data = mean(data(:,freq_min_idx:freq_max_idx), 2); % avg over the selected freq range
        
        % append the data from all subjects in this group
        data_for_anova = vertcat(data_for_anova, data); 
        % also append this group's label accordingly (same number of times as how many subjects were appended)
        for count = 1:length(data)
            grouping_var = [grouping_var; groups{g}]; 
        end
    end

    [p,tbl,stats] = anova1(data_for_anova, grouping_var, 'off');
    F = cell2mat(tbl(2,5)); % get the F-value
    % p-value is already returned; it's the same as cell2mat(tbl(2,6))
    
    % print results to console output
    disp([freq_band ' band:']);
    disp(['F = ', num2str(F)]) % print the F-value
    disp(['p = ', num2str(p)]) % print the p-value
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== Individual channel analysis ===== %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ANOVA at each channel & each freq (27 x 30 = 810 comparisons)
% (see Figure 2a in Flavia paper)

load([stats_folder groups{1} '\GA_individuals.mat']); % load one group as an example
N_chan = size(GA_freq_indi.powspctrm, 2); % number of channels
N_freq = size(GA_freq_indi.powspctrm, 3); % number of freqs

% initialise "chan x freq" matrix to store f-values & p-values
f_values = zeros(N_chan, N_freq);
p_values = zeros(N_chan, N_freq);

for i = 1:N_chan % loop through each channel
    for j = 1:N_freq % loop through each freq
        
        % collate data for anova (each subject should have a single value)
        % https://au.mathworks.com/help/stats/anova1.html
        % https://au.mathworks.com/matlabcentral/answers/159736-how-do-i-perform-unbalanced-anova
        data_for_anova = [];
        grouping_var = {};
        for g = 1:length(groups)
            load([stats_folder groups{g} '\GA_individuals.mat']);
            if logged % apply log transformation if needed
                GA_freq_indi.powspctrm = log(GA_freq_indi.powspctrm);
            end
            % each subject has a single value, representing the power at this chan & this freq
            data = GA_freq_indi.powspctrm(:,i,j);
            
            % append the data from all subjects in this group
            data_for_anova = vertcat(data_for_anova, data); 
            % also append this group's label accordingly (same number of times as how many subjects were appended)
            for count = 1:length(data)
                grouping_var = [grouping_var; groups{g}]; 
            end
        end
        
        [p,tbl,stats] = anova1(data_for_anova, grouping_var, 'off');
        
        % store the F-value & p-value into the "chan x freq" matrix
        f_values(i,j) = cell2mat(tbl(2,5));
        p_values(i,j) = cell2mat(tbl(2,6));
    end
end

save([stats_folder_indi_chan 'stats--chan_x_freq.mat'], 'f_values', 'p_values');

%% plots
% using a separate section here as the ANOVAs above take a while to run,
% so we just load the stats results previously saved
load([stats_folder_indi_chan 'stats--chan_x_freq.mat']);

figure; title('F-values');
imagesc(f_values)
colorbar
ylabel('EEG channel');
xlabel('Frequency (Hz)');

export_fig(gcf, [stats_folder_indi_chan 'indi-chan-analysis' logged_suffix '_F-values.png']);

% not necessary to plot the p-values (it would just be the opposite of the f-values)
%{
figure; title('p-values');
imagesc(p_values, [0 0.05]) % only plot p-values up to 0.05
colorbar
ylabel('EEG channel');
xlabel('Frequency (Hz)');

export_fig(gcf, [stats_folder_indi_chan 'indi-chan-analysis' logged_suffix '_p-values.png']);
%}


%% ANOVA at each channel for each freq band (27 x 3 = 81 comparisons)
% (see Figure 2b in Flavia paper)

% loop thru each freq band
for band = 1:length(FREQ_BANDS)
    freq_band = cell2mat(FREQ_BANDS{band, 1}); % first field is the freq band name
    freq_range = FREQ_BANDS{band, 2}; % second field is the freq range in Hz

    % find the start index & end index for the selected freq range
    load([stats_folder groups{1} '\GA_individuals.mat']); % load one group as an example
    freq_min_idx = find(GA_freq_indi.freq == freq_range(1));
    freq_max_idx = find(GA_freq_indi.freq == freq_range(end));
    
    
    % initialise an array to store the f-value for each channel
    N_chan = size(GA_freq_indi.powspctrm, 2); % get the number of channels
    f_values = zeros(N_chan, 1);
    p_values = zeros(N_chan, 1);

    for i = 1:N_chan % loop through each channel
        % collate data for anova (each subject should have a single value - see above)
        data_for_anova = [];
        grouping_var = {};
        for g = 1:length(groups)
            load([stats_folder groups{g} '\GA_individuals.mat']);
            if logged % apply log transformation if needed
                GA_freq_indi.powspctrm = log(GA_freq_indi.powspctrm);
            end
            data = GA_freq_indi.powspctrm(:,i,freq_min_idx:freq_max_idx); % extract the power values at this channel (only for the selected freq range)
            data = mean(data,3); % take the avg over that freq range (each subject now ends up with a single value)
            
            % append the data from all subjects in this group
            data_for_anova = vertcat(data_for_anova, data); 
            % also append this group's label accordingly (same number of times as how many subjects were appended)
            for count = 1:length(data)
                grouping_var = [grouping_var; groups{g}]; 
            end
        end

        [p,tbl,stats] = anova1(data_for_anova, grouping_var, 'off');
        f_values(i) = cell2mat(tbl(2,5));
        p_values(i) = cell2mat(tbl(2,6));
    end

    % create a dummy var for plotting
    freq = GA_freq;
    freq.powspctrm = f_values;
    freq.freq = 0; % we no longer have a frequency dimension, just fill with a dummy value

    % plot topography based on the f-values
    plot_TFR_topo(freq, lay, freq_band, [], [stats_folder_indi_chan 'fvalues_' logged_suffix '_']);
end


%% Cluster-based statistical analysis
% https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/
% https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_stats/#2-compute-between-participants-contrasts

% NOTE: At the moment this only works if you are comparing 3 groups!
%       If you have more groups, the code below needs to be manually adjusted.


load(NEIGHBOURS_FILE); % obtained using 'triangulation' method in ft_prepare_neighbour

% load the data
load([stats_folder groups{1} '\allSubjects_freq.mat']);
group1 = allSubjects_freq;
load([stats_folder groups{2} '\allSubjects_freq.mat']);
group2 = allSubjects_freq;
load([stats_folder groups{3} '\allSubjects_freq.mat']);
group3 = allSubjects_freq;

if logged % using logged power
    for i = 1:length(group1)
        group1{i}.powspctrm = log(group1{i}.powspctrm); % apply log transformation
    end
    for i = 1:length(group2)
        group2{i}.powspctrm = log(group2{i}.powspctrm);
    end
    for i = 1:length(group3)
        group3{i}.powspctrm = log(group3{i}.powspctrm);
    end
end

% For future - if we want to try automatically adjusting for number of groups
%{
for g = 1:length(groups)
    load([stats_folder groups{g} '\allSubjects_freq.mat']);
    if logged % using logged power
        for i = 1:length(allSubjects_freq)
            allSubjects_freq{i}.powspctrm = log(allSubjects_freq{i}.powspctrm); % apply log transformation
        end
    end
end
%}


% Opt 1: find spatio-freq cluster:
%freq_band = 'spatio-freq';
%freq_range = [1 30];
% Opt 2: find spatial cluster for a particular freq band:
freq_band = 'alpha';
freq_range = [9 12];

cfg = [];
cfg.channel = 'all';
cfg.frequency = freq_range;
if strcmp(freq_band, 'spatio-freq')
    cfg.avgoverfreq = 'no'; % if searching for spatio-freq cluster, don't avg
else
    cfg.avgoverfreq = 'yes'; % if using a particular freq band, avg over that band
end

cfg.method = 'montecarlo';
cfg.correctm = 'cluster'; %'no'; % it is common in MEG studies to run uncorrected at cfg.alpha = 0.001
cfg.clusteralpha = 0.05; % threshold for selecting candidate samples to form clusters
cfg.clusterstatistic = 'maxsum';
%cfg.clusterthreshold = 'nonparametric_common';
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment
cfg.minnbchan = 2; % minimum number of neighbourhood channels required to be significant 
                   % in order to form a cluster 
                   % (default: 0, ie. each single channel can be considered a cluster)

cfg.alpha = 0.05; %0.001  % threshold for cluster-level statistics (any cluster with a p-value lower than this will be reported as sig - an entry of '1' in .mask field)
cfg.numrandomization = 1000; % Rule of thumb: use 500, and double this number if it turns out 
    % that the p-value differs from the chosen alpha (e.g. 0.05) by less than 0.02

% design matrix
num_group1 = length(group1);
num_group2 = length(group2);
num_group3 = length(group3);

design = zeros(1, num_group1 + num_group2 + num_group3);
design(1, 1:num_group1) = 1;
design(1, (num_group1 + 1):(num_group1 + num_group2)) = 2;
design(1, (num_group1 + num_group2 + 1):end) = 3;

cfg.design = design;
cfg.ivar   = 1;

cfg.statistic = 'indepsamplesF'; % independent samples F-test (3 groups of subjects)
%cfg.tail = 1;
cfg.clustertail = 1;

[stat] = ft_freqstatistics(cfg, group1{:}, group2{:}, group3{:});
length(find(stat.mask)) % display how many chans were significant/marginal

save([stats_folder_cluster 'minnbchan' mat2str(cfg.minnbchan) '_' freq_band '.mat'], 'stat');


%% ft_clusterplot
% this is a wrapper around ft_topoplot, automatically extracts info about the cluster

stats_filename = [stats_folder_cluster 'minnbchan2_' freq_band];
load([stats_filename '.mat']);

% use a nice-looking colourmap
%ft_hastoolbox('brewermap', 1); % ensure this toolbox is on the path
%cmap = colormap(flipud(brewermap(64, 'RdBu')));

% too much warning, can't see the console output
ft_warning off 'FieldTrip:ft_clusterplot:ft_topoplotTFR:topoplot_common:ft_selectdata:getdimord:warning_dimord_could_not_be_determined:line681'

cfg = [];
%cfg.zlim = [-5 5]; % set scaling (range of t-values) (usually using automatic is ok) 
cfg.highlightcolorpos = [0 0 0]; % black for pos clusters
cfg.highlightcolorneg = [255/255 192/255 203/255]; % pink for neg clusters
cfg.alpha = 0.05; % any clusters with a p-value below this threshold will be plotted
cfg.layout = lay;
cfg.style = 'straight';
%cfg.colormap = cmap;

% turn on the following lines if you are after one particular subplot
cfg.subplotsize = [1 1];
cfg.highlightsizeseries = [9 9 9 9 9]; % make the highlight markers bigger
cfg.colorbar = 'yes'; % shows the scaling

ft_clusterplot(cfg, stat);

export_fig(gcf, [stats_filename '.png']);
