path = 'Z:\PRJ-Transient\MIGRAINES\EEG\Subject 676\Session 1\';
filename = '20150716172032.edf';

% define trials
cfg            = [];
cfg.dataset    = [path filename];
cfg.continuous = 'yes';
cfg.channel    = 'all';
data           = ft_preprocessing(cfg);

% visually inspect the data
cfg           = [];
cfg.viewmode  = 'vertical';
cfg.continous = 'yes';
cfg.blocksize = 120; % display 2-min segments
cfg.ylim      = [ -32   32 ];
ft_databrowser(cfg, data);