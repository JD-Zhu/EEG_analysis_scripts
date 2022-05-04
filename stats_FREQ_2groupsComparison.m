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
% (note - 2 groups is a special case, should use t-test rather than F-test;
% but can try & see if the anova1 fn can actually handle 2 groups)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global LAYOUT_FILE; global NEIGHBOURS_FILE; global PLOT_XLIM; global is_conn;
global FREQ_BANDS;
common();

load(LAYOUT_FILE);


% = Settings =

% Run t-tests using equal variance (i.e. pooled) or unequal variance? (p-values obtained were similar)
% 'unequal' is the most conservative approach
% https://www.investopedia.com/terms/t/t-test.asp
% https://statmagic.info/Content/Help-Content/two-sample-mean.html
varType = 'unequal'; %'equal';

% Use logged power or absolute power?
logged = true;

% PLEASE SPECIFY the folder for this statistical analysis
stats_folder = 'Z:\Analysis\Macefield Lab\Bec\stats\4vs4\';


% where to read in the freq results for each group:
migraineurs_folder = [stats_folder '\GA_migraineurs\'];
controls_folder = [stats_folder '\GA_controls\'];

% load data
load([migraineurs_folder 'GA_avg.mat']);
mig_avg = GA_freq;
load([controls_folder 'GA_avg.mat']);
ctrl_avg = GA_freq;

load([migraineurs_folder 'GA_individuals.mat']);
mig_indi = GA_freq_indi; 
load([controls_folder 'GA_individuals.mat']);
ctrl_indi = GA_freq_indi;


%% for standard freq analysis
if ~is_conn

    % apply log transformation if needed & set appropriate file paths
    if logged
        mig_indi.powspctrm_abs = mig_indi.powspctrm; % retain a copy of the absolute power (just in case)
        mig_indi.powspctrm = log(mig_indi.powspctrm); % apply log transformation
        ctrl_indi.powspctrm_abs = ctrl_indi.powspctrm; % retain a copy of the absolute power (just in case)
        ctrl_indi.powspctrm = log(ctrl_indi.powspctrm); % apply log transformation

        logged_suffix = '_logged';
    else
        logged_suffix = '';
    end


    % locations to save results
    stats_folder_indi_chan = [stats_folder '\indi-chan-analysis\'];
    mkdir(stats_folder_indi_chan);
    stats_folder_cluster = [stats_folder_indi_chan '\cluster_stat' logged_suffix '\'];
    mkdir(stats_folder_cluster);


    %% plot overall power (mig vs ctrl)
    
    % absolute power
    figure; hold on;
    plot(mig_avg.freq, mean(mig_avg.powspctrm));
    plot(ctrl_avg.freq, mean(ctrl_avg.powspctrm));

    xlim(PLOT_XLIM);
    xlabel('Frequency (Hz)');
    ylabel('Absolute power (uV^2)');
    legend({'Migraineurs', 'Controls'});
    hold off;

    export_fig(gcf, [stats_folder 'overall_power_mig-vs-ctrl.png']);

    % log-transformed power
    figure; hold on;
    plot(mig_avg.freq, mean(log(mig_avg.powspctrm)));
    plot(ctrl_avg.freq, mean(log(ctrl_avg.powspctrm)));

    xlim(PLOT_XLIM);
    xlabel('Frequency (Hz)');
    ylabel('Power (log[uV^2]');
    legend({'Migraineurs', 'Controls'});
    hold off;

    % add shaded region for alpha band
    %{
    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
    x = [9 12 12 9]; % specify x,y coordinates of the 4 corners
    y = [ylow ylow yhigh yhigh]; 
    alpha = 0.3; % use alpha to set transparency
    patch(x,y,'black', 'FaceAlpha',alpha, 'EdgeColor','none', ...
        'HandleVisibility','off') % turn off HandleVisibility so it won't show up in the legends
    ylim(ylimits); % ensure ylim doesn't get expanded
    %}
    
    export_fig(gcf, [stats_folder 'overall_power_logged_mig-vs-ctrl.png']);


    %% Analysis of overall power
    % run t-test (mig vs ctrl) for each freq band
    
    % loop thru each freq band
    for band = 1:length(FREQ_BANDS)
        freq_band = cell2mat(FREQ_BANDS{band, 1}); % first field is the freq band name
        freq_range = FREQ_BANDS{band, 2}; % second field is the freq range in Hz
    
        % find the start index & end index for the selected freq range
        freq_min_idx = find(ctrl_indi.freq == freq_range(1));
        freq_max_idx = find(ctrl_indi.freq == freq_range(end));
        
        
        % avg over all channels
        mig_overall = squeeze(mean(mig_indi.powspctrm, 2)); 
        ctrl_overall = squeeze(mean(ctrl_indi.powspctrm, 2));
        
        % avg over the selected freq range
        a = mean(mig_overall(:,freq_min_idx:freq_max_idx), 2); 
        b = mean(ctrl_overall(:,freq_min_idx:freq_max_idx), 2);

        
        % this function conducts a two-sample t-test
        [h,p,ci,stats] = ttest2(a, b, 'Vartype',varType); % or should it be equal variance? (p-values were similar)
        
        % print results to console output
        disp([freq_band ' band:']);
        stats.tstat  % t-value
        p            % p-value
    end


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ===== Individual channel analysis ===== %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% t-test at each channel & each freq (27 x 30 = 810 comparisons)
    % (see Figure 2a in Flavia paper)

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
    imagesc(t_values, [-4.5 4.5]) % adjust clims manually
    colorbar
    ylabel('EEG channel');
    xlabel('Frequency (Hz)');

    export_fig(gcf, [stats_folder_indi_chan 'indi-chan-analysis' logged_suffix '.png']);


    %% t-test at each channel for each freq band (27 x 3 = 81 comparisons)
    % (see Figure 2b in Flavia paper)

    % loop thru each freq band
    for band = 1:length(FREQ_BANDS)
        freq_band = cell2mat(FREQ_BANDS{band, 1}); % first field is the freq band name
        freq_range = FREQ_BANDS{band, 2}; % second field is the freq range in Hz

        % find the start index & end index for the selected freq range
        freq_min_idx = find(ctrl_indi.freq == freq_range(1));
        freq_max_idx = find(ctrl_indi.freq == freq_range(end));
        
    
        N_chan = size(mig_indi.powspctrm, 2); % get the number of channels

        % initialise an array to store the t-value for each channel
        t_values = zeros(N_chan, 1);
        p_values = zeros(N_chan, 1);

        for i = 1:N_chan % loop through each channel
            a = mig_indi.powspctrm(:,i,freq_min_idx:freq_max_idx); % extract power for all migraineurs (only for the selected freq range)
            a = mean(a,3); % take the avg over that freq range
            b = ctrl_indi.powspctrm(:,i,freq_min_idx:freq_max_idx); % do the same for controls
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
        zlim = [0 3];
        plot_TFR_topo(freq, lay, freq_band, [], [stats_folder_indi_chan 'tvalues_' logged_suffix '_'], zlim);
    end
    

    %% plot the topography difference btwn two groups
    % e.g. subtract the topo for alpha (migraineurs minus controls)

    % NOTE: this plot is based on absolute power diff in the GA, 
    % whereas Flavia plotted the t-values (see above - "Figure 2b")

    % calculate the difference GA
    diff_GA = mig_avg;
    diff_GA.powspctrm = mig_avg.powspctrm - ctrl_avg.powspctrm;

    % plot topography based on the difference GA
    for band = 1:length(FREQ_BANDS)
        freq_band = cell2mat(FREQ_BANDS{band, 1}); % first field is the freq band name
        freq_range = FREQ_BANDS{band, 2}; % second field is the freq range in Hz
        plot_TFR_topo(diff_GA, lay, freq_band, [freq_range(1) freq_range(end)], [stats_folder 'topo_mig-minus-ctrl_'])
    end
    

    %% Cluster-based statistical analysis
    % https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_freq/
    % https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_stats/#2-compute-between-participants-contrasts

    % load the data
    load([migraineurs_folder 'allSubjects_freq.mat']);
    mig = allSubjects_freq;
    load([controls_folder 'allSubjects_freq.mat']);
    ctrl = allSubjects_freq;

    load(NEIGHBOURS_FILE); % obtained using 'triangulation' method in ft_prepare_neighbour

    if logged % using logged power
        for i = 1:length(mig)
            mig{i}.powspctrm_abs = mig{i}.powspctrm; % retain a copy of the absolute power (just in case)
            mig{i}.powspctrm = log(mig{i}.powspctrm); % apply log transformation
        end
        for i = 1:length(ctrl)
            ctrl{i}.powspctrm_abs = ctrl{i}.powspctrm;
            ctrl{i}.powspctrm = log(ctrl{i}.powspctrm);
        end
    end


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

    cfg.alpha = 0.1; %0.001  % threshold for cluster-level statistics (any cluster with a p-value lower than this will be reported as sig - an entry of '1' in .mask field)
    cfg.numrandomization = 1000; % Rule of thumb: use 500, and double this number if it turns out 
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

    save([stats_folder_cluster 'minnbchan' mat2str(cfg.minnbchan) '_' freq_band '.mat'], 'stat');


    %% ft_clusterplot (plots t-values by default)
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
    cfg.alpha = 0.1; % any clusters with a p-value below this threshold will be plotted
    cfg.layout = lay;
    cfg.style = 'straight';
    %cfg.colormap = cmap;

    % turn on the following lines if you are after one particular subplot
    cfg.subplotsize = [1 1];
    cfg.highlightsizeseries = [9 9 9 9 9]; % make the highlight markers bigger
    cfg.colorbar = 'yes'; % shows the scaling

    ft_clusterplot(cfg, stat);

    export_fig(gcf, [stats_filename '.png']);

else

    %% %%%%%%%%%%%%%%%%%%%%%%%%%
    % ===== Connectivity ===== %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    The cohspctrm is a 4D matrix, e.g. for mig_indi, it is
    12x27x27x30  (subj_chan_chan_freq)

    Qs: 
    1. does it make sense to avg the coh values across a freq band, like what 
    I did below? (there are no other alternatives atm)
    L: yes.

    2. the coh values are percentages / scores between 0 ~ 1, is it appropriate
    to run t-tests on such data?
    L: yes if the values are normally distributed.
    %}

    % for some reason mig #2 data is messed up, so remove for now
    %mig_indi.cohspctrm(2,:,:,:) = [];

    
    %% Compare coherence btwn mig & ctrl (t-test for each pair of channels: 27 x 27)
    % the actual number of tests (i.e. non-repeated channel pairs) is 1 + 2 + ... + 26 = 351 comparisons

    % loop thru each freq band
    for band = 1:length(FREQ_BANDS)
        freq_band = cell2mat(FREQ_BANDS{band, 1}); % first field is the freq band name
        freq_range = FREQ_BANDS{band, 2}; % second field is the freq range in Hz

        % find the start index & end index for the selected freq range
        freq_min_idx = find(ctrl_indi.freq == freq_range(1));
        freq_max_idx = find(ctrl_indi.freq == freq_range(end));

    
        N_chan = size(mig_indi.cohspctrm, 2); % number of channels

        c = []; % for normality check (see below)

        % initialise "chan x chan" matrix to store t-values
        t_values = zeros(N_chan, N_chan);
        p_values = zeros(N_chan, N_chan);

        for i = 1:N_chan % loop through each "from" channel
            for j = 1:N_chan % loop through each "to" channel
                if i > j % only test each channel pair once (i.e. triangular matrix, not square)
                    a = mig_indi.cohspctrm(:,i,j, freq_min_idx:freq_max_idx); % extract coh values for all migraineurs
                    a = mean(a, 4); % avg over the selected freq range
                    b = ctrl_indi.cohspctrm(:,i,j, freq_min_idx:freq_max_idx); % extract coh values for all controls
                    b = mean(b, 4); % avg over the selected freq range

                    % check if data are normally distributed (an assumption of the t-test)
                    % use Anderson-Darling test for normality: 0 = normally distributed; 1 = not normally distributed
                    c = [c adtest([a; b])]; % opt 1: test normality on the data for each channel pair, and collect results into an array
                    %c = [c; a; b]; % opt 2: put all data together, you can then run adtest(c) afterwards
                    
                    [h,p,ci,stats] = ttest2(a, b, 'Vartype',varType);
                    t_values(i,j) = stats.tstat; % store the t-value into the "chan x chan" matrix


                    % if not normally distributed, try non-parametric alternatives of t-test:
                    % https://www.statisticshowto.com/probability-and-statistics/non-normal-distributions/

                    % (1) Wilcoxon rank sum test / Mann-Whitney U-test
                    [p,h,stats] = ranksum(a, b, 'tail','right'); % right-tailed test: a > b
                    % (2) Kruskal-Wallis test (replacement for one-way ANOVA)
                    %grouping_var = [repmat({'pt'}, [length(a),1]); repmat({'ctrl'}, [length(b),1])];
                    %[p,tbl,stats] = kruskalwallis([a; b], grouping_var, 'off');

                    p_values(i,j) = p; % store the p-value into the "chan x chan" matrix
                    
                else % no test done in these cells, mark them as NaN (so they can be plotted as white colour below)
                    t_values(i,j) = NaN;
                    p_values(i,j) = NaN;
                end
            end
        end

        length(find(c)) % out of 351 possible channel pairs, how many have non-normally distributed data
        %adtest(c)

        
        %t_values(isnan(t_values)) = -1; % replace all NaNs first, so they don't show up as sig in the plot
        figure; title('t-values (migraineurs > controls)');
        imagesc(t_values, 'AlphaData', ~isnan(t_values), [-3 3]); % any cells that were marked as NaN (i.e. no test done) will be plotted with a transparency alpha of 0 (i.e. completely transparent)
                                                                  % set the clims manually so that all plots are of same scaling
        colorbar
        ylabel('EEG channel');
        xlabel('EEG channel');
        export_fig(gcf, [stats_folder 'tvalues_' freq_band '.png']);

        %p_values(isnan(p_values)) = 1; % replace all NaNs first, so they don't show up as sig in the plot
        figure;
        imagesc(p_values, 'AlphaData', ~isnan(p_values), [0 0.05]) % only plot p-values up to 0.05
        colorbar
        ylabel('EEG channel');
        xlabel('EEG channel');
        export_fig(gcf, [stats_folder 'pvalues_MannWhitneyUtest_migGreater_' freq_band '.png']);
        %export_fig(gcf, [stats_folder 'pvalues_kruskalwallis_' freq_band '.png']);
    end
end