path = 'Z:\PRJ-Transient\MIGRAINES\EEG\controls\Subject_891\Session 1 (32 Channel)\';
      %'Z:\GRP-Henderson\RawPOWMRI_MRIData01\SCITrigemMRIData\Subject_700\Session 1 (32 Channel)\EEG\';
file = dir([path '*.edf']);

% define trials
cfg            = [];
cfg.dataset    = [path file.name];
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

% visually inspect the data
cfg           = [];
cfg.viewmode  = 'vertical';
cfg.continous = 'yes';
cfg.blocksize = 120; % display 2-min segments
cfg.ylim      = [ -1600  1600 ];
ft_databrowser(cfg, data);