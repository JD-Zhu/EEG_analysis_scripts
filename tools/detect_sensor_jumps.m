%
% Author: Judy Zhu (github.com/JD-Zhu)
%

% select according to which system you are using: adult/child
ref_channels = 161:163;

% specify confile location
rawfile = 'E:\Judy\Exp2\6_MEG-data\RAW_DATA\B02-YW-3523\3523_ME180_2019_08_29_B1-concat_TSPCA10000_3.con';    


% Read in the data
cfg                      = [];
cfg.trialfun             = 'ft_trialfun_general';
cfg.headerfile           = rawfile;
cfg.datafile             = rawfile;
cfg.trialdef.triallength = Inf;
cfg.trialdef.ntrials     = 1;
cfg = ft_definetrial(cfg);

cfg.continuous = 'yes';


%% Option 1: use ft_artifact_jump (a wrapper around ft_artifact_zvalue, 
% with tailored settings for detecting 'jump' artefact)

cfg.artfctdef.jump.channel = ref_channels; %'all'; %'MEG';
cfg.artfctdef.jump.cutoff = 8;
cfg.artfctdef.jump.interactive = 'yes';

[cfg, artifact_jump] = ft_artifact_jump(cfg);


%% Option 2: use ft_artifact_zvalue (choose your own settings)
%{
% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel = ref_channels; %'all'; %'MEG';
cfg.artfctdef.zvalue.cutoff = 8;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.cumulative = 'yes';
cfg.artfctdef.zvalue.medianfilter = 'yes';
cfg.artfctdef.zvalue.medianfiltord = 9;
cfg.artfctdef.zvalue.absdiff = 'yes';

% make the process interactive
cfg.artfctdef.zvalue.interactive = 'yes';

[cfg, artifact_jump] = ft_artifact_zvalue(cfg);
%}