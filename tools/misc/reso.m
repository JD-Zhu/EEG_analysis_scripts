%
% Author: Judy Zhu (github.com/JD-Zhu)
%

% select the raw file here
rawfile = 'E:\Judy\Exp2_pilots\Pilot_stage2\Data\resting-state\HPF0.3_gain5.con';            

cfg                      = [];
cfg.trialfun             = 'ft_trialfun_general';
cfg.headerfile           = rawfile;
cfg.datafile             = rawfile;
cfg.trialdef.triallength = Inf;
cfg.trialdef.ntrials     = 1; % read in all data as a single segment, coz filtering should be done on continuous data
cfg = ft_definetrial(cfg);

cfg.continuous = 'yes';
alldata = ft_preprocessing(cfg);



% now search through all the measured values
M = alldata.trial{1,1};

% initialise the var, will be used to store the reso for each channel
reso = zeros(1, 160);

% loop thru each channel
for i = 1:160
    X = M(i,:); % extract all values for this channel

    X = sort(X);
    step = diff(X); % find the diff btwn each 2 adjacent values
    step(find(step == 0)) = []; % remove the 0s (same measured values)
    
    reso(i) = min(step); % find the min change from one measurement to the next
end

% output the minimum step we found
min(reso)