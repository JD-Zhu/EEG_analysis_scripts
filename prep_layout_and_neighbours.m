% Only need to run this if we are using a new MEG system,
% or the sensor layout changes


% Create layout file for later (plotting)
cfg      = [];
cfg.grad = alldata.grad; % struct containing gradiometer definition
lay      = ft_prepare_layout(cfg, alldata); % creates a 2-D layout of the channel locations
save([ResultsFolder 'lay.mat'], 'lay');


% Prepare neighbours & save for use in stat analysis
cfg_neighb        = [];
cfg_neighb.method = 'triangulation'; % 'distance' may have some issues
           % 'triangulation' seems to give each channel more neighbours 
           % (about twice as many as 'distance' method gives)
neighbours        = ft_prepare_neighbours(cfg_neighb, all_blocks);
save([ResultsFolder 'neighbours.mat'], 'neighbours');