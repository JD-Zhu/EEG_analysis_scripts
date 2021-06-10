%{
Plot topography showing location of the sensors which formed the cluster.

ft_clusterplot is convenient for early visualisation of significant channel 
locations. It is a wrapper around ft_topoplot, and automatically extracts
the spatial and temporal extent of the cluster from the stat output. It 
produces a series of subplots, spanning the time window in which the effect
was found (1 subplot per sample).

If you do not wish to include all the subplots in your paper, you can produce 
a plot averaged over time points (i.e. this script). You need to make the 
decision on how to include the spatial extent (i.e. electrodes) of your
cluster. You could show all electrodes that were part of the cluster at any
one time point, or you could plot the electrode size as a function of how
often (i.e. at how many time points) the electrode was part of the cluster.

Here we plot the topography averaged over the interval of the effect (425-550ms), 
and highlight the channels which are significant at the "middle" time sample (485ms).

See here for more information:
http://www.fieldtriptoolbox.org/faq/how_not_to_interpret_results_from_a_cluster-based_permutation_test
%}

%% run the #define section to obtain values for global vars
global ResultsFolder; 
common();

% SELECT which set of single-subject ERFs to use
run_name = 'TSPCA10000_3'; % this should be a folder name inside the "Results_ERF" folder
ResultsFolder_thisrun = [ResultsFolder run_name '\\STATS_axial###\\']; % ERF results for all subjects

% load the results
load([ResultsFolder_thisrun 'stats_Interactions_minnbchan4.mat']);
load('lay.mat');

% SELECT which effect to plot & SPECIFY the effect duration
stat_output = MixCost_interaction; % 155-200ms (minnbchan = 2); 160-180ms (minnbchan = 3); 165-170ms (minnbchan = 4);
start_time = 0.175;
end_time = 0.175;

% load nice colourmap
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64, 'RdBu')));


%% define what contrast you want to plot & compute the GA difference
% (only need to run this section if plotting topography based on GA, i.e. Opt 1 below)

load([ResultsFolder_thisrun 'GA_avg.mat']); % takes a long time to load

% first, define the 2 conds to be compared
% choose one of these:

% (a) if you want 3 separate plots (ie. pairwise comparisons, showing mixing effect in Nat, Art & Bi)
% then modify the following accordingly for the 3 contexts (each time you 
% need to run the code below again, up to the part where the topoplot is generated)
GA_1 = GA_erf.NatStay;
GA_2 = GA_erf.NatSingle;

%{
% (b) if you want to plot the mixing-effect diff btwn 2 contexts (ie. interaction)
% e.g. here we plot MixCost_nat_vs_bi
GA_1 = GA_erf.BiStay;
GA_1.avg = GA_erf.BiStay.avg - GA_erf.BiSingle.avg; % mixing cost in Bi
GA_2 = GA_erf.NatStay;
GA_2.avg = GA_erf.NatStay.avg - GA_erf.NatSingle.avg; % mixing cost in Nat
%}

% then, calc the diff btwn the 2 conds
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
timelock = ft_math(cfg, GA_1, GA_2);


%% plot topography for axial gradiometers

cfg                   = [];
cfg.layout            = lay;
cfg.colormap          = cmap;
cfg.colorbar          = 'yes'; % shows the scales
%cfg.colorbar         = 'EastOutside';
%cfg.zlim              = [0 10]; %'maxabs'; % set the scaling

cfg.xlim = [start_time end_time]; % duration of the effect (as reported by ft_clusterplot)
                          % topography will be averaged over this interval

cfg.highlight       = 'on'; % highlight significant channels
cfg.highlightsymbol = '.';
cfg.highlightsize   = 26;
cfg.highlightcolor  = 'w';


% which channels to highlight? 

% here we highlight all channels that were sig 
% at one or more time points during the effect interval
sig_channels = [];
% each cycle checks one channel
for i = 1:size(stat_output.mask, 1)
    if find(stat_output.mask(i,:)) % this channel was sig at some point
        sig_channels = [sig_channels i]; % so add it to the list
    end
end
cfg.highlightchannel  = stat_output.label(sig_channels);

% Alternatively, we highlight the channels that are sig at the middle time point
%{
time_point = (start_time + end_time) / 2; 
time_point = time_point * 1000; % convert unit to ms
%time_point = 440; % set manually if needed
time_index = round(time_point / 5 + 1); % index of the time point you want, round it up/down to an integer
                                        % note: this is only correct if the epoch starts from time 0;
                                        % if it starts from -100ms, then need to adjust index accordingly
%}
% set manually if needed
%time_index = 56 or 57 (minnbchan = 2); 55 (minnbchan = 3); 
%cfg.highlightchannel  = stat_output.label(ismember(stat_output.posclusterslabelmat(:,time_index),1)); % find all the channels showing '1' at that time point


cfg.style             = 'straight'; % no contour lines, only the colour gradients
cfg.comment           = 'no';
cfg.gridscale         = 512;
cfg.marker            = 'on'; % show the location of all channels?


% Opt 1: plot topography based on actual erf amplitude (i.e. GA)
%cfg.parameter = 'avg'; 
%ft_topoplotER(cfg, timelock);

% Opt 2: plot topography based on F-values in the stat output
cfg.parameter = 'stat'; % what field to plot, e.g. selecting 'stat' will plot t-values / F-values
ft_topoplotER(cfg, stat_output);

set(gca,'fontsize',22); % set colorbar fontsize
%title(gca, 'F-values');


%% transform into planar gradient

% Opt 1: planar transformation on actual erf amplitude

%TODO% apply planar transformation on the GA computed above.
%Q: is it valid to transform the GA diff? or need to transform each GA, then calc the diff again?
% or even need to do transformation on each subject's ERF first, then compute GA and then GA diff?


% Opt 2: planar transformation on F-values in stat output <-- is this even valid??
%{
% make a fake GA_erf structure, then plug in the F-values from the stat output
timelock = GA_erf.NatStay; % copy the GA structure from anywhere
timelock.avg = stat_output.stat; % not sure if it uses the 'avg' or 'var' field to calc planar,
timelock.var = stat_output.stat; % so I'll just replace both
timelock.time = stat_output.time;
%}


% fill in the "grad" field so that ft_megplanar will work
timelock.grad = stat_output.grad; 

% Compute the planar gradient at each sensor location
% in both the horizontal and the vertical direction 
cfg                 = [];
cfg.feedback        = 'yes';
cfg.method          = 'template';
load('neighbours_tri.mat');
cfg.neighbours      = neighbours;

cfg.planarmethod    = 'sincos';
planar              = ft_megplanar(cfg, timelock); % plug in any timelock struct

% Combine the horizontal and vertical components 
% of the planar gradient using Pythagoras' Rule
planarComb = ft_combineplanar([], planar);

% Plot the planar gradient
cfg = [];
%cfg.zlim = [0 3];
cfg.xlim = [start_time end_time]; % duration of the effect (as reported by ft_clusterplot)
                          % topography will be averaged over this interval
cfg.layout = lay;
cfg.colormap = cmap;
cfg.colorbar = 'yes';

cfg.style             = 'straight'; % no contour lines, only the colour gradients
cfg.feedback          = 'no'; % necessary to specify this when cfg.style = 'straight', otherwise errors!!
cfg.comment           = 'no';
cfg.gridscale         = 512;
cfg.marker            = 'off'; % show the location of all channels?

ft_topoplotER(cfg, planarComb);
set(gca,'fontsize',22); % set colorbar fontsize
