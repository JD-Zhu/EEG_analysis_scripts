% This utility is useful when you see something strange in the ERF (e.g.
% unexpected peak before 0, strange drifts, etc.)
% Use this script to browse trial by trial, to see if it was a few 
% crazy epochs that caused the average to look crazy
%
% Author: Judy Zhu (github.com/JD-Zhu)
%

%% select the trials you want to check
cond_trials = trials_clean.NatStay;


%% browse trial by trial
cfg           = [];
cfg.viewmode  = 'butterfly';
cfg.ylim      = [ -4e-13   4e-13 ];

ft_databrowser(cfg, cond_trials);


%% remove bad trials

% manually enter the indices of bad trials below
bad_trials = [41 38 50 51 61 33 34 35 16 11 10 9 8 4 3 2 1];

good_trials = setdiff(1:length(cond_trials.trial), bad_trials);
cfg             = [];
cfg.trials      = good_trials;
cond_trials_clean = ft_redefinetrial(cfg, cond_trials);


%% compute & plot ERF
erf = ft_timelockanalysis([], cond_trials_clean);

cfg              = [];
cfg.showlabels   = 'yes';
cfg.fontsize     = 6;
cfg.layout       = lay;
cfg.baseline     = ERF_BASELINE; % makes no diff if we've already done baseline correction earlier
cfg.baselinetype = 'absolute';
cfg.graphcolor   = cell2mat(colours); 
cfg.linestyle    = lineTypes;
cfg.xlim = PLOT_XLIM;

figure;
ft_multiplotER(cfg, erf);
    