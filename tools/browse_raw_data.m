%
% Author: Judy Zhu (github.com/JD-Zhu)
%

% select the raw file here
rawfile = 'E:\Judy\Exp2\6_MEG-data\RAW_DATA\A04-ZC-3545\3545_ME180_2019_09_04_B1-concat_TSPCA10000_3.con';            

cfg                      = [];
cfg.trialfun             = 'ft_trialfun_general';
cfg.headerfile           = rawfile;
cfg.datafile             = rawfile;
cfg.trialdef.triallength = Inf;
cfg.trialdef.ntrials     = 1; % read in all data as a single segment, coz filtering should be done on continuous data
cfg = ft_definetrial(cfg);

cfg.continuous = 'yes';
alldata = ft_preprocessing(cfg);

% select MEG channels
%{
cfg         = [];
cfg.channel = alldata.label(1:160);
%cfg.channel = alldata.label([1:160 191 193]);
%cfg.channel = alldata.label([113 117 145 146 148 150 153 155 156]);
alldata = ft_selectdata(cfg, alldata);
%}

% Try detrend
%{
cfg = []; % demean == remove 0-order polynomial
cfg.detrend = 'yes'; % == remove 1st-order polynomial (i.e. linear trend)
cfg.polyremoval = 'yes'; % == remove higher order polynomial (default: 2nd-order)
cfg.polyorder = 22; % customise the order (higher order => affects higher freqs?)
alldata_detrend = ft_preprocessing(cfg, alldata);
%}

% Plot
cfg           = [];
cfg.viewmode  = 'vertical';
cfg.continous = 'yes';
cfg.blocksize = 60; % display 60-sec segments
cfg.ylim      = [ -4e-13   4e-13 ];
% randomly sample a few channels, to check if the bandpass filter 
% can successfully filter out the low-freq drift
%{
%cfg.channel = [1 8 15 23 30 37 44 51 58 65 72 78 84 91 98 104 111 119 125 132 139 146 153]; 
cfg.channel = [1 11 35 49 60 72 85 98 111 125 139 153]; 
cfg.ylim = [ -1.0941e-12  1.0941e-12 ];
%}
ft_databrowser(cfg, alldata);

% Plot after detrend
%ft_databrowser(cfg, alldata_detrend);
