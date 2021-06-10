% After running stats_ROI.m, you can use this script to perform post-hoc
% t-tests btwn 2 conds in a particular ROI.
%
% Author: Judy Zhu (github.com/JD-Zhu)
%

% SELECT THE ROI (where you saw the interaction effect):
ROI_name = 'LIFG';
fprintf(['\nROI: ' ROI_name '\n']);

data = allSubjects_ROIs.(ROI_name); % data for the current ROI

% set some config for the statistical test
cfg = [];
cfg.channel = {'all'}; % there is only one channel (i.e. the virtual sensor for this ROI)
cfg.avgoverchan = 'yes'; % this is necessary (or else FT will ask for cfg.neighbours)

cfg.latency = [-0.1 0.75]; % time interval over which the experimental 
                     % conditions must be compared (in seconds)
cfg.avgovertime = 'no'; % if yes, this will average over the entire time window chosen in cfg.latency 
                        % (useful when you want to look at a particular component, e.g. to look at M100,
                        % cfg.latency = [0.08 0.12]; cfg.avgovertime = 'yes'; )

%load([ResultsFolder_ROI 'neighbours.mat']); % this is the sensor layout - it's the same for all subjects (even same across experiments). So just prepare once & save, then load here
%cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT'; %cfg.statistic = 'ft_statfun_indepsamplesT'; OR 'ft_statfun_depsamplesFmultivariate';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
%cfg.minnbchan = 3; % minimum number of neighbourhood channels required to be significant 
                   % in order to form a cluster 
                   % (default: 0, ie. each single channel can be considered a cluster).
                   % 4 or 5 is a good choice; 2 is too few coz it's even below
                   % the resolution of the sensor layout(??)

cfg.tail = 0;
cfg.clustertail = 0; % 2 tailed test
cfg.alpha = 0.1; % report all effects with p < 0.1
cfg.correcttail = 'prob'; % correct for 2-tailedness
cfg.numrandomization = 2000; % Rule of thumb: use 500, and double this number if it turns out 
    % that the p-value differs from the critical alpha-level (0.05 or 0.01) by less than 0.02

numSubjects = length(files);
within_design_2x2 = zeros(2, 2*numSubjects);
within_design_2x2(1, :) = repmat(1:numSubjects, 1, 2);
within_design_2x2(2, 1:numSubjects) = 1;
within_design_2x2(2, numSubjects+1:2*numSubjects) = 2;

cfg.design = within_design_2x2;
cfg.uvar  = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar  = 2; % row of design matrix that contains independent variable (i.e. the conditions)


% Run the t-tests

fprintf('\nCUE window -> Testing switch cost in English:\n');
[stat_Eng] = ft_timelockstatistics(cfg, data.cueenstay{:}, data.cueenswitch{:});
fprintf('\nCUE window -> Testing switch cost in Chinese:\n');
[stat_Chn] = ft_timelockstatistics(cfg, data.cuechstay{:}, data.cuechswitch{:});

fprintf('\nCUE window -> Stay trials: English vs Chinese:\n');
[stat_stay] = ft_timelockstatistics(cfg, data.cueenstay{:}, data.cuechstay{:});
fprintf('\nCUE window -> Switch trials: English vs Chinese:\n');
[stat_switch] = ft_timelockstatistics(cfg, data.cueenswitch{:}, data.cuechswitch{:});
