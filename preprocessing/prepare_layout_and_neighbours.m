% only need to do this once & save for later


%% prepare layout

% read the sensor positions
% https://www.fieldtriptoolbox.org/faq/how_are_electrodes_magnetometers_or_gradiometers_described/
%chanpos = readtable('chanlocs_XYZ_61_nolabel.txt');
%elec.chanpos = chanpos;
%elec.label = all_labels;
%elec.unit = 'mm';
load('elec.mat'); % just load the version we have already made, this contains 62 channels (with ref sensor CPz added back)

cfg = [];
cfg.elec = elec;
lay = ft_prepare_layout(cfg);
figure; ft_plot_layout(lay);

save lay lay


%% prepare neighbours
%{
Default distance is 4.
For MQ system:
  default distance gives on average 4.6 neighbours per channel;       (prob too few)
  distance of 5 gives on average 8.0 neighbours per channel;
  distance of 6 gives on average 11.4 neighbours per channel;
  distance of 7 gives on average 15.0 neighbours per channel;         
  'triangulation' method gives on average 7.6 neighbours per channel. (I used this for MEG Exp 1)
%}

cfg        = [];
cfg.layout = lay;
cfg.method = 'triangulation'; %'distance'; 
%cfg.neighbourdist = 7; % distance threshold (only applicable to 'distance' method)
cfg.feedback = 'yes'; % this will produce a plot showing the neighbour relationships
% http://www.fieldtriptoolbox.org/faq/how_can_i_define_neighbouring_sensors/

neighbours = ft_prepare_neighbours(cfg); 
save neighbours neighbours
