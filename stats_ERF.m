%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_ERF.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Grand average & statistical analysis on sensor-space ERFs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
%clear all;


% = Settings =
% Please adjust as required:

% SELECT which set of single-subject ERFs to use
run_name = 'TSPCA10000_3'; % this should be a folder name inside the "Results_ERF" folder

% Apply planar transformation before running stats?
% (this will make everything +ve, therefore avoiding the sign-flipping issue)
% But it's not necessarily more valid than the original (ie. axial) version,
% as the more processing you do to the data, the more chance for distortion.
PLANAR_TRANSFORM = false;

% perform channel repair on each subject's ERF?
CHANNEL_REPAIR = false; % only need to do the repair once, we'll save repaired erf in the folder below
repaired_erf_folder = 'channelrepaired\\'; % need to create this folder first

% cfg.avgovertime setting in cluster-based permutation test
AVGOVERTIME = false;
%TIME_WINDOW_TO_AVG = [0.060 0.110]; % must set this var, if AVGOVERTIME is set to true


%%
% run the #define section
global ResultsFolder; % all subjects' erf data are stored here
global filename_suffix; % erf results file suffix

global eventnames_real; global colours_and_lineTypes; 
global colours; global lineTypes;
global PLOT_XLIM; global ERF_BASELINE;
common();

% location of ERF results for all subjects
if (PLANAR_TRANSFORM)
    ResultsFolder_thisrun = [ResultsFolder run_name '\\STATS_planar\\'];
else
    ResultsFolder_thisrun = [ResultsFolder run_name '\\STATS_axial\\'];
end

% initialise allSubjects_erf (each field holds all subjects' erf in that condition)
allSubjects_erf.NatStay = {};
allSubjects_erf.NatSwitch = {};
allSubjects_erf.NatSingle = {};
allSubjects_erf.ArtStay = {};
allSubjects_erf.ArtSwitch = {};
allSubjects_erf.ArtSingle = {};
allSubjects_erf.BiStay = {};
allSubjects_erf.BiSwitch = {};
allSubjects_erf.BiSingle = {};


%% Read data

% SKIP THIS SECTION - The output have now been saved ("allSubjects_erf.mat" in "STATS_axial" folder).
% Data will be loaded in relevant sections below.

% find all .mat files in ResultsFolder_thisrun
files = dir([ResultsFolder_thisrun '*_erf' filename_suffix '.mat']);

% each cycle reads in one '.mat' file (ie. one subject's erf data)
for i = 1:length(files)
    filename = [ResultsFolder_thisrun files(i).name];
    load(filename);
        
    for j = 1:length(eventnames_real) % 9 conds (if collapsed across langs) or 18 conds (if not collapsed)
        % perform channel repair if needed
        if (CHANNEL_REPAIR == true)
            load('neighbours.mat');
            load('all_labels.mat');
            erf_clean.(eventnames_real{j}) = repair_bad_channels(erf_clean.(eventnames_real{j}), neighbours, all_labels);
        end
        % add to allsubjects matrix
        allSubjects_erf.(eventnames_real{j}) = [allSubjects_erf.(eventnames_real{j}) erf_clean.(eventnames_real{j})];
    end
    
    % save the new erf after channel repair
    if (CHANNEL_REPAIR == true)
        save([ResultsFolder_thisrun repaired_erf_folder files(i).name], 'SubjectFolder', 'erf_clean');
    end
end


%% Descriptives
% http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock#within-subjects_experiments

% SKIP THIS SECTION - The output have now been saved ("GA_avg.mat" & "GA_individuals.mat").
% Data will be loaded in relevant sections below.

fprintf('\n= COMPUTING & PLOTTING CROSS-SUBJECT AVERAGES =\n');

% CALCULATE the grand average (across all subjects) for each condition
cfg = [];
cfg.channel   = {'all', '-AG101', '-AG122', '-AG007', '-AG103'}; % remove noisy sensors (MEG Exp2):
                                                                % ch100 (AG101) is always noisy -> Remove for all ptps!
                                                                % ch006 & ch102 also shows the same noise occasionally.
                                                                % ch121 (AG122) tends to show square noise.
cfg.latency   = 'all';
cfg.parameter = 'avg';
for j = 1:length(eventnames_real)
    cfg.keepindividual = 'no'; % average across subjects
    GA_erf.(eventnames_real{j}) = ft_timelockgrandaverage(cfg, allSubjects_erf.(eventnames_real{j}){:});  

    cfg.keepindividual = 'yes'; % do not average across subjects, keep the data for each individual subject
    GA_indi.(eventnames_real{j}) = ft_timelockgrandaverage(cfg, allSubjects_erf.(eventnames_real{j}){:}); 

    % "{:}" means to use data from all elements of the variable
end

% save the GA files
GA_output_file = [ResultsFolder_thisrun 'GA_avg.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_erf');
end
GA_output_file = [ResultsFolder_thisrun 'GA_individuals.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_indi');
end

% multiplot
load('lay.mat');
        
cfg              = [];
cfg.showlabels   = 'yes';
cfg.fontsize     = 6;
cfg.layout       = lay;
cfg.baseline     = ERF_BASELINE;
cfg.baselinetype = 'absolute';
cfg.graphcolor   = cell2mat(colours); 
cfg.linestyle    = lineTypes;
cfg.xlim         = PLOT_XLIM;

figure('Name','ft_multiplotER: GA_erf (9 conds)');
% convert struct to cell array, then you can feed it in as 'varargin'
cellarray = struct2cell(GA_erf);
ft_multiplotER(cfg, cellarray{:});
% specify the legends manually (otherwise it will display incorrectly)
lines = findall(gcf, 'Type','line');
lines = lines([29 26 23 20 17 14 11 8 5]); % grab the correct lines (this is complicated because many lines are plotted in ft_multiplot)
legend(lines, eventnames_real);


% CALCULATE global averages across all sensors (i.e. GFP = global field power)
cfg        = [];
cfg.method = 'power';
%cfg.channel = {'AG017', 'AG018', 'AG019', 'AG022', 'AG023', 'AG025', 'AG029', 'AG063', 'AG064', 'AG143'}; % 10 sig channels in cluster
cfg.channel   = {'all', '-AG101', '-AG122', '-AG007', '-AG103'}; % remove noisy sensors (see above)
for j = 1:length(eventnames_real)
    GA_erf_GFP.(eventnames_real{j}) = ft_globalmeanfield(cfg, GA_erf.(eventnames_real{j}));
end

% plot GFP
figure('Name','GFP_all_subjects'); hold on
for j = 1:length(eventnames_real)
    if colours_and_lineTypes % use a combination of colours and line types to distinguish conds
        plot(GA_erf_GFP.(eventnames_real{j}).time, GA_erf_GFP.(eventnames_real{j}).avg, 'color', colours{j}, 'LineStyle', lineTypes{j});
    else % just use colours
        plot(GA_erf_GFP.(eventnames_real{j}).time, GA_erf_GFP.(eventnames_real{j}).avg, 'color', colours(j,:));
    end
    xlim(PLOT_XLIM);
end
legend(eventnames_real);


% average across all 4 conds (for selecting windows for peaks)
averageAcrossConds = GA_erf_GFP.NatStay;
averageAcrossConds.avg = (GA_erf_GFP.NatStay.avg + GA_erf_GFP.NatSwitch.avg + GA_erf_GFP.NatSingle.avg ...
                        + GA_erf_GFP.ArtStay.avg + GA_erf_GFP.ArtSwitch.avg + GA_erf_GFP.ArtSingle.avg ...
                        + GA_erf_GFP.BiStay.avg + GA_erf_GFP.BiSwitch.avg + GA_erf_GFP.BiSingle.avg) / 9;

figure('Name','GFP_all_subjects - Averaged across all conds'); 
plot(averageAcrossConds.time, averageAcrossConds.avg); 
xlim([-0.2 0.8]);


%% Apply planar transformation before running stats?
% (this will make everything +ve, therefore avoiding the sign-flipping issue)

% SKIP THIS SECTION - The output have now been saved ("allSubjects_erf.mat" in "STATS_planar" folder).
% Data will be loaded in relevant sections below.

if (PLANAR_TRANSFORM)
    
    load('neighbours.mat');
    
    counter = 0; % just to check the correct # of ERFs are processed

    cond_names = fieldnames(allSubjects_erf);
    for j = 1:length(cond_names) % each cycle handles one cond (e.g. NatStay)
        cond_name = cond_names{j};

        for i = 1:length(allSubjects_erf.(cond_name)) % each cycle handles one subject
            % Compute the planar gradient at each sensor location
            % in both the horizontal and the vertical direction 
            cfg                 = [];
            cfg.feedback        = 'no';
            cfg.method          = 'template';
            cfg.planarmethod    = 'sincos';
            cfg.neighbours      = neighbours;

            planar = ft_megplanar(cfg, allSubjects_erf.(cond_name){i});

            % Combine the horizontal and vertical components 
            % of the planar gradient using Pythagoras' Rule
            allSubjects_erf.(cond_name){i} = ft_combineplanar([], planar);

            counter = counter + 1;
        end
    end
    
    fprintf(['\nPlanar transformed ' int2str(counter) ' ERFs (should be 9*24=216).\n']);
end


%% Statistical analysis

% load the data
load([ResultsFolder_thisrun 'allSubjects_erf.mat']); 

data = allSubjects_erf; % change to an easy name
clear allSubjects_erf; % clear up the memory


% Optional: select which subjects to (not) use
for j = 1:length(eventnames_real)
    %data.(eventnames_real{j})([1 11 16 19 20 21 14]) = [];
end


%%
fprintf('\n= STATS: CLUSTER-BASED PERMUTATION TESTS =\n');

cfg = [];
cfg.channel   = {'all', '-AG101', '-AG122', '-AG007', '-AG103'}; % remove noisy sensors (see above)
load('neighbours_tri.mat'); % obtained using 'trigangulation' method in ft_prepare_neighbour
%load('neighbours_dist6.mat'); % obtained using 'distance' method, with a distance threshold of 6 (default is 4)
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

% can choose diff time windows to analyse for cue epochs & target epochs
% (these will be fed into cfg.latency accordingly)
if (AVGOVERTIME)
    latency_cue = TIME_WINDOW_TO_AVG; % time range to average over
    cfg.avgovertime = 'yes'; % if yes, this will average over the entire time window chosen in cfg.latency 
                            % (useful when you want to look at a particular component, e.g. to look at M100,
                            % cfg.latency = [0.08 0.12]; cfg.avgovertime = 'yes'; )
else % autoly detect temporal cluster
    latency_cue = [-0.1 0.6]; % time interval over which the experimental 
                               % conditions must be compared (in seconds)
    cfg.avgovertime = 'no';
end

cfg.method = 'montecarlo';
cfg.correctm = 'cluster'; %'no'; % it is common in MEG studies to run uncorrected at cfg.alpha = 0.001
cfg.clusteralpha = 0.05; % threshold for selecting candidate samples to form clusters
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2; % minimum number of neighbourhood channels required to be significant 
                   % in order to form a cluster 
                   % (default: 0, ie. each single channel can be considered a cluster).
                   % 4 or 5 is a good choice; 2 is too few coz it's even below
                   % the resolution of the sensor layout (i.e. 2 adjacent sensors might
                   % really be measuring the same thing, so ofc they are both sig)

cfg.alpha = 0.05; %0.001  % threshold for cluster-level statistics (any cluster with a p-value lower than this will be reported as sig - an entry of '1' in .mask field)
cfg.numrandomization = 2000; % Rule of thumb: use 500, and double this number if it turns out 
    % that the p-value differs from the chosen alpha (e.g. 0.05) by less than 0.02

numSubjects = length(data.(eventnames_real{1})); % check how many subjects we are including
within_design_1x2 = zeros(2,2*numSubjects);
within_design_1x2(1,:) = repmat(1:numSubjects,1,2);
within_design_1x2(2,1:numSubjects) = 1;
within_design_1x2(2,numSubjects+1:2*numSubjects) = 2;

within_design_1x3 = zeros(2, 3*numSubjects);
within_design_1x3(1, :) = repmat(1:numSubjects, 1, 3);
within_design_1x3(2, 1:numSubjects) = 1;
within_design_1x3(2, numSubjects+1:2*numSubjects) = 2;
within_design_1x3(2, 2*numSubjects+1:3*numSubjects) = 3;

cfg.uvar  = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar  = 2; % row of design matrix that contains independent variable (i.e. the conditions)

cfg.latency = latency_cue;

% 4 time windows to avg over
%{
% peak latencies: 90, 150, 245, 495
cfg.latency = [0.060 0.110];
cfg.latency = [0.130 0.175];
cfg.latency = [0.200 0.280]; % bi sw$ (corresponding to the RdlPFC effect in ### version) can be found in [0.180 0.265] or [0.170 0.290];
cfg.latency = [0.350 0.640]; % bi sw$ is found here (for cfg.minnbchan = 2)
%}


%----- Run the statistical tests -----%

%% Use 'F test' for interaction & main effects (coz there are 3 levels in "context")
cfg.statistic = 'ft_statfun_depsamplesFunivariate';
cfg.design = within_design_1x3;
cfg.tail = 1; % -1 = left, 1 = right
cfg.clustertail = 1; % for F test, can only select right-sided tail
                     % https://github.com/fieldtrip/fieldtrip/blob/master/statfun/ft_statfun_depsamplesFunivariate.m

% INTERACTIONS
cfg.minnbchan = 2; % from my experience so far, generally 2~4 finds the
                   % same cluster (2 has the widest temporal span & largest
                   % number of channels in the cluster, while 4 is on the 
                   % other end). At 5, no clusters are found.

% sw$ interaction (i.e. calc sw$ in each context, then submit the 3 sw$ to F-test)
fprintf('\n= Sw$ interaction (i.e. compare sw$ in 3 contexts using an F test) =\n');
[timelock_SwCost_Nat, timelock_SwCost_Bi] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatStay, data.NatSwitch, data.BiStay, data.BiSwitch); %'2-1 vs 4-3'
[timelock_SwCost_Art, ~] = combine_conds_for_T_test('fieldtrip', 'interaction', data.ArtStay, data.ArtSwitch, data.BiStay, data.BiSwitch);
[SwCost_interaction] = ft_timelockstatistics(cfg, timelock_SwCost_Nat{:}, timelock_SwCost_Art{:}, timelock_SwCost_Bi{:});

% mix$ interaction (i.e. calc mix$ in each context, then submit the 3 mix$ to F-test)
fprintf('\n= Mix$ interaction (i.e. compare mix$ in 3 contexts using an F test) =\n');
[timelock_MixCost_Nat, timelock_MixCost_Bi] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatSingle, data.NatStay, data.BiSingle, data.BiStay);
[timelock_MixCost_Art, ~] = combine_conds_for_T_test('fieldtrip', 'interaction', data.ArtSingle, data.ArtStay, data.BiSingle, data.BiStay);
[MixCost_interaction] = ft_timelockstatistics(cfg, timelock_MixCost_Nat{:}, timelock_MixCost_Art{:}, timelock_MixCost_Bi{:});

length(find(SwCost_interaction.mask))
length(find(MixCost_interaction.mask))

%save([ResultsFolder_thisrun 'stats_Interactions_minnbchan' mat2str(cfg.minnbchan) '.mat'], 'SwCost_interaction', 'MixCost_interaction');


% MAIN EFFECT of context (collapsed across single-stay-switch)
fprintf('\n= Main effect of context (F test) =\n');
timelock_Nat = data.NatStay;
timelock_Art = data.ArtStay;
timelock_Bi = data.BiStay;
for i = 1:numSubjects
    timelock_Nat{i}.avg = (data.NatStay{i}.avg + data.NatSwitch{i}.avg + data.NatSingle{i}.avg) / 3; % average across stay/switch/single
    timelock_Art{i}.avg = (data.ArtStay{i}.avg + data.ArtSwitch{i}.avg + data.ArtSingle{i}.avg) / 3; % average across stay/switch/single
    timelock_Bi{i}.avg = (data.BiStay{i}.avg + data.BiSwitch{i}.avg + data.BiSingle{i}.avg) / 3; % average across stay/switch/single
end
[Main_Context] = ft_timelockstatistics(cfg, timelock_Nat{:}, timelock_Art{:}, timelock_Bi{:});
    
length(find(Main_Context.mask))
    
% MAIN EFFECT of ttype (collapsed across nat-art-bi)
fprintf('\n= Main effect of ttype (F test) =\n');
timelock_Stay = data.NatStay;
timelock_Switch = data.NatSwitch;
timelock_Single = data.NatSingle;
for i = 1:numSubjects
    timelock_Stay{i}.avg = (data.NatStay{i}.avg + data.ArtStay{i}.avg + data.BiStay{i}.avg) / 3; % average across nat/art/bi
    timelock_Switch{i}.avg = (data.NatSwitch{i}.avg + data.ArtSwitch{i}.avg + data.BiSwitch{i}.avg) / 3; 
    timelock_Single{i}.avg = (data.NatSingle{i}.avg + data.ArtSingle{i}.avg + data.BiSingle{i}.avg) / 3;
end
[Main_Ttype] = ft_timelockstatistics(cfg, timelock_Single{:}, timelock_Stay{:}, timelock_Switch{:});
    
length(find(Main_Ttype.mask))

% Alternatively,
% we can test "switch effect" & "mixing effect" separately (no need to correct for 2 comparisons)
cfg.statistic = 'depsamplesT'; % t-test (i.e. for comparing 2 conds)
cfg.design = within_design_1x2;
% make sure we are using 2-tailed for t-test (using 1-tailed is generally frowned upon)
cfg.tail = 0; % -1 = left, 1 = right, 0 = 2-tailed
cfg.clustertail = 0; 
cfg.correcttail = 'prob'; % correct for 2-tailedness

fprintf('\nMain effect of switch: (t-test)\n');
[Switch] = ft_timelockstatistics(cfg, timelock_Stay{:}, timelock_Switch{:});
fprintf('\nMain effect of mix: (t-test)\n');
[Mix] = ft_timelockstatistics(cfg, timelock_Single{:}, timelock_Stay{:});

length(find(Switch.mask))             % not sig
length(find(Mix.mask))                % sig (did not survive Bonferroni at minnbchan = 3; survived at minnbchan = 2)

%save([ResultsFolder_thisrun 'stats_MainEffects_minnbchan' mat2str(cfg.minnbchan) '.mat'], 'Main_Context', 'Main_Ttype', 'Switch', 'Mix');


%% PLANNED PAIRWISE COMPARISONS within each context 
% (previously known as "SANITY CHECK")
fprintf('\n\n= Planned pairwise comparisons to assess sw$ & mix$ within each context\n');
    
% Make sure we are using 2-tailed t-tests:
cfg.statistic = 'depsamplesT'; % t-test (i.e. for comparing 2 conds)
cfg.design = within_design_1x2;
cfg.tail = 0;
cfg.clustertail = 0; % 2 tailed test
cfg.correcttail = 'prob'; % correct for 2-tailedness

% Switch cost in each context
[Nat_sw] = ft_timelockstatistics(cfg, data.NatStay{:}, data.NatSwitch{:}); %allSubj_cue_ch_switchCost{:}, allSubj_cue_en_switchCost{:});
[Art_sw] = ft_timelockstatistics(cfg, data.ArtStay{:}, data.ArtSwitch{:});
[Bi_sw] = ft_timelockstatistics(cfg, data.BiStay{:}, data.BiSwitch{:}); 

% Mixing cost in each context
[Nat_mix] = ft_timelockstatistics(cfg, data.NatSingle{:}, data.NatStay{:});
[Art_mix] = ft_timelockstatistics(cfg, data.ArtSingle{:}, data.ArtStay{:});
[Bi_mix] = ft_timelockstatistics(cfg, data.BiSingle{:}, data.BiStay{:}); 

length(find(Nat_sw.mask))  % not sig
length(find(Art_sw.mask))  % sig (but did not survive Bonferroni)
length(find(Bi_sw.mask))   % not sig
length(find(Nat_mix.mask)) % sig (but did not survive Bonferroni)
length(find(Art_mix.mask)) % not sig
length(find(Bi_mix.mask))  % sig (but did not survive Bonferroni)

%save([ResultsFolder_thisrun 'stats_pairwise_minnbchan' mat2str(cfg.minnbchan) '.mat'], 'Nat_sw', 'Art_sw', 'Bi_sw', 'Nat_mix', 'Art_mix', 'Bi_mix');


%% UNPACKING main effects & interactions

cfg.statistic = 'depsamplesT'; % t-test (i.e. for comparing 2 conds)
cfg.design = within_design_1x2;
% make sure we are using 2-tailed for t-test (using 1-tailed is generally frowned upon)
cfg.tail = 0; % -1 = left, 1 = right, 0 = 2-tailed
cfg.clustertail = 0; 
cfg.correcttail = 'prob'; % correct for 2-tailedness


% Unpack mix$ interaction
stat = MixCost_interaction;

% read out the sensors in the cluster
[row,col] = find(stat.mask); 
sig_chans = unique(row)';

% Take plain avg over all sig sensors
cfg.channel = stat.label(sig_chans);
cfg.avgoverchan = 'yes'; 
cfg.minnbchan = 0;      % if avgoverchan = 'yes', then set this to 0, otherwise you cannot 
                        % possibly form a cluster (because there is only one "channel")
% Also take avg over time
cfg.latency = [0.155 0.200]; % duration of the cluster
cfg.avgovertime = 'yes';


% run relevant code (should have been run during F-test above) to compute these:
% timelock_MixCost_Nat, timelock_MixCost_Art, timelock_MixCost_Bi
[timelock_MixCost_Nat, timelock_MixCost_Bi] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatSingle, data.NatStay, data.BiSingle, data.BiStay);
[timelock_MixCost_Art, ~] = combine_conds_for_T_test('fieldtrip', 'interaction', data.ArtSingle, data.ArtStay, data.BiSingle, data.BiStay);

% test the 2x2 contrasts
[MixCost_nat_vs_bi] = ft_timelockstatistics(cfg, timelock_MixCost_Nat{:}, timelock_MixCost_Bi{:});
[MixCost_art_vs_bi] = ft_timelockstatistics(cfg, timelock_MixCost_Art{:}, timelock_MixCost_Bi{:});
[MixCost_nat_vs_art] = ft_timelockstatistics(cfg, timelock_MixCost_Nat{:}, timelock_MixCost_Art{:});

length(find(MixCost_nat_vs_bi.mask))  % sig
length(find(MixCost_art_vs_bi.mask))  % not sig
length(find(MixCost_nat_vs_art.mask)) % sig

% test the pairwise contrast (single vs stay) in each context
[Nat_mix_posthoc] = ft_timelockstatistics(cfg, data.NatSingle{:}, data.NatStay{:});
[Art_mix_posthoc] = ft_timelockstatistics(cfg, data.ArtSingle{:}, data.ArtStay{:});
[Bi_mix_posthoc] = ft_timelockstatistics(cfg, data.BiSingle{:}, data.BiStay{:}); 

length(find(Nat_mix_posthoc.mask)) % sig
length(find(Art_mix_posthoc.mask)) % not sig
length(find(Bi_mix_posthoc.mask))  % sig

%save([ResultsFolder_thisrun 'stats_Interactions_minnbchan' mat2str(cfg.minnbchan) '_unpack_avgovertime.mat'], 'MixCost_nat_vs_bi', 'MixCost_art_vs_bi', 'MixCost_nat_vs_art', 'Nat_mix_posthoc', 'Art_mix_posthoc', 'Bi_mix_posthoc');

% We now use avgovertime/avgoverchan to unpack main effects & interactions,
% so the code below is obsolete.
%{
% TO UNPACK THE INTERACTIONS, we compare the sw$ & mix$ for each pair of contexts (i.e. 3 t-tests)
% unpack sw$ interaction
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatStay, data.NatSwitch, data.BiStay, data.BiSwitch);
[SwCost_nat_vs_bi] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.ArtStay, data.ArtSwitch, data.BiStay, data.BiSwitch);
[SwCost_art_vs_bi] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatStay, data.NatSwitch, data.ArtStay, data.ArtSwitch);
[SwCost_nat_vs_art] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});

length(find(SwCost_nat_vs_bi.mask))  % not sig
length(find(SwCost_art_vs_bi.mask))  % not sig
length(find(SwCost_nat_vs_art.mask)) % not sig

% unpack mix$ interaction
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatSingle, data.NatStay, data.BiSingle, data.BiStay);
[MixCost_nat_vs_bi] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); 
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.ArtSingle, data.ArtStay, data.BiSingle, data.BiStay);
[MixCost_art_vs_bi] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatSingle, data.NatStay, data.ArtSingle, data.ArtStay);
[MixCost_nat_vs_art] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); 

length(find(MixCost_nat_vs_bi.mask)) % sig (but did not survive Bonferroni)
length(find(MixCost_art_vs_bi.mask)) % not sig
length(find(MixCost_nat_vs_art.mask)) % not sig

save([ResultsFolder_thisrun 'stats_Interactions_unpack_minnbchan' mat2str(cfg.minnbchan) '.mat'], 'SwCost_nat_vs_bi', 'SwCost_art_vs_bi', 'SwCost_nat_vs_art', 'MixCost_nat_vs_bi', 'MixCost_art_vs_bi', 'MixCost_nat_vs_art');


% TO UNPACK THE MAIN EFFECT of context
fprintf('\nMain effect of context - unpacking (3 t-tests):');
timelock_Nat = data.NatStay;
timelock_Art = data.ArtStay;
timelock_Bi = data.BiStay;
for i = 1:numSubjects
    timelock_Nat{i}.avg = (data.NatStay{i}.avg + data.NatSwitch{i}.avg + data.NatSingle{i}.avg) / 3; % average across stay/switch/single
    timelock_Art{i}.avg = (data.ArtStay{i}.avg + data.ArtSwitch{i}.avg + data.ArtSingle{i}.avg) / 3; % average across stay/switch/single
    timelock_Bi{i}.avg = (data.BiStay{i}.avg + data.BiSwitch{i}.avg + data.BiSingle{i}.avg) / 3; % average across stay/switch/single
end
%cfg.latency = latency_cue; % time interval over which the experimental 
fprintf('\n  -> Nat vs Bi');
[Context_nat_vs_bi] = ft_timelockstatistics(cfg, timelock_Nat{:}, timelock_Bi{:});
fprintf('\n  -> Art vs Bi');
[Context_art_vs_bi] = ft_timelockstatistics(cfg, timelock_Art{:}, timelock_Bi{:});
fprintf('\n  -> Nat vs Art');
[Context_nat_vs_art] = ft_timelockstatistics(cfg, timelock_Nat{:}, timelock_Art{:});

length(find(Context_nat_vs_bi.mask))  % sig (but did not survive Bonferroni)
length(find(Context_art_vs_bi.mask))  % not sig
length(find(Context_nat_vs_art.mask)) % not sig

save([ResultsFolder_thisrun 'stats_MainEffects_unpack_minnbchan' mat2str(cfg.minnbchan) '.mat'], 'Context_nat_vs_bi', 'Context_art_vs_bi', 'Context_nat_vs_art', 'Switch', 'Mix');
%}


%% Below are from MEG Exp 1
% Interaction (i.e. calc sw$ in each lang, then test the 2 sw$)
% http://www.fieldtriptoolbox.org/faq/how_can_i_test_an_interaction_effect_using_cluster-based_permutation_tests
%{
% manual calculation above is now replaced by combine_conds_for_T_Test()
fprintf('\nCUE window -> Testing lang x ttype interaction:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
cfg.latency = latency_cue;  
[cue_interaction] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_cue_ch_switchCost{:}, allSubj_cue_en_switchCost{:});
fprintf('\nTARGET window -> Testing lang x ttype interaction:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
cfg.latency = latency_target;  
[target_interaction] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_target_ch_switchCost{:}, allSubj_target_en_switchCost{:});

% Main effect of lang (collapse across stay-switch)
%{
for i = 1:numSubjects
    allSubj_cue_ch{i} = allSubjects_erf.cuechstay{i};
    allSubj_cue_ch{i}.avg = (allSubjects_erf.cuechstay{i}.avg + allSubjects_erf.cuechswitch{i}.avg) / 2; % cue_ch_all.avg = (cuechstay.avg + cuechsw.avg) / 2
    allSubj_cue_en{i} = allSubjects_erf.cueenstay{i};
    allSubj_cue_en{i}.avg = (allSubjects_erf.cueenstay{i}.avg + allSubjects_erf.cueenswitch{i}.avg) / 2;
    allSubj_target_ch{i} = allSubjects_erf.targetchstay{i};
    allSubj_target_ch{i}.avg = (allSubjects_erf.targetchstay{i}.avg + allSubjects_erf.targetchswitch{i}.avg) / 2; % target_ch_all.avg = (targetchstay.avg + targetchsw.avg) / 2
    allSubj_target_en{i} = allSubjects_erf.targetenstay{i};
    allSubj_target_en{i}.avg = (allSubjects_erf.targetenstay{i}.avg + allSubjects_erf.targetenswitch{i}.avg) / 2;    
end
%}
fprintf('\nCUE window -> Main effect of lang:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_12vs34', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
cfg.latency = latency_cue; 
[cue_lang] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_cue_ch{:}, allSubj_cue_en{:});
fprintf('\nTARGET window -> Main effect of lang:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_12vs34', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
cfg.latency = latency_target;
[target_lang] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_target_ch{:}, allSubj_target_en{:});

% Main effect of switch (collapse across langs)
%{
for i = 1:numSubjects
    allSubj_cue_stay{i} = allSubjects_erf.cuechstay{i};
    allSubj_cue_stay{i}.avg = (allSubjects_erf.cuechstay{i}.avg + allSubjects_erf.cueenstay{i}.avg) / 2;
    allSubj_cue_switch{i} = allSubjects_erf.cuechswitch{i};
    allSubj_cue_switch{i}.avg = (allSubjects_erf.cuechswitch{i}.avg + allSubjects_erf.cueenswitch{i}.avg) / 2;
    allSubj_target_stay{i} = allSubjects_erf.targetchstay{i};
    allSubj_target_stay{i}.avg = (allSubjects_erf.targetchstay{i}.avg + allSubjects_erf.targetenstay{i}.avg) / 2;
    allSubj_target_switch{i} = allSubjects_erf.targetchswitch{i};
    allSubj_target_switch{i}.avg = (allSubjects_erf.targetchswitch{i}.avg + allSubjects_erf.targetenswitch{i}.avg) / 2;
end
%}
fprintf('\nCUE window -> Main effect of ttype:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_13vs24', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
cfg.latency = latency_cue; 
[cue_ttype] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_cue_stay{:}, allSubj_cue_switch{:});
fprintf('\nTARGET window -> Main effect of ttype:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_13vs24', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
cfg.latency = latency_target; 
[target_ttype] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_target_stay{:}, allSubj_target_switch{:});

% check for effects by searching the .mask field
effect_cue_interaction = length(find(cue_interaction.mask)) % if not 0, then we have an effect here
effect_cue_lang = length(find(cue_lang.mask))
effect_cue_ttype = length(find(cue_ttype.mask))
effect_target_interaction = length(find(target_interaction.mask))
effect_target_lang = length(find(target_lang.mask))
effect_target_ttype = length(find(target_ttype.mask))

%save([ResultsFolder_thisrun 'stats.mat'], 'cue_interaction', 'cue_lang', 'cue_ttype', 'target_interaction', 'target_lang', 'target_ttype');
%}


%% Plotting: use ft_clusterplot & ft_topoplot

%load([ResultsFolder_thisrun 'stats.mat']);
load('lay.mat');

% use a nice-looking colourmap
ft_hastoolbox('brewermap', 1); % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64, 'RdBu')));

% select which comparison to plot
stat = MixCost_interaction; % here we plot the only effect that seems to survive correction (at minnbchan = 0)
                  % to explore where (both in terms of time & location) the effect might have possibly
                  % occurred
                  % [TODO] then we can define more precise time window &
                  % set avgovertime = 'yes', which should give us more
                  % sensitivity, and allow us to increase the minnbchan to
                  % a reasonable number: 2 (ft tutorial) or 4 (Paul)

%% ft_clusterplot (plots t-values by default)
% this is a wrapper around ft_topoplot, automatically extracts info about the cluster

% too much warning, can't see the console output
ft_warning off 'FieldTrip:ft_clusterplot:ft_topoplotTFR:topoplot_common:ft_selectdata:getdimord:warning_dimord_could_not_be_determined:line681'

cfg = [];
%cfg.zlim = [-5 5]; % set scaling (range of t-values) (usually using automatic is ok) 
cfg.highlightcolorpos = [1 1 1]; % white for pos clusters
cfg.highlightcolorneg = [255/255 192/255 203/255]; % pink for neg clusters
cfg.alpha = 0.05; % any clusters with a p-value below this threshold will be plotted
cfg.layout = lay;
cfg.colormap = cmap;

% turn on the following lines if you are after one particular subplot
%cfg.subplotsize = [1 1];
%cfg.colorbar = 'yes'; % shows the scaling

ft_clusterplot(cfg, stat);

% for some reason, the plots still use old colormap
% running the following again sometimes updates the plots to new colormap
ft_hastoolbox('brewermap', 1); % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64, 'RdBu')));


%% ft_topoplot (can plot t-values or actual erf amplitude) 
% need to write everything yourself (ie. pretty much re-implement whats's in ft_clusterplot).
% if u don't want to do this manually, might be able to set cfg.parameter = GA.avg when calling ft_clusterplot,
% to make use of the automatic processing provided by ft_clusterplot
%{
load([ResultsFolder_thisrun 'GA_avg.mat']);

% first, define the 2 conds to be compared (this time using cross-subject averages, i.e. GA)
% here we look at main effect of ttype in cue window, so we collapse across langs
GA_cue_stay = GA_erf.cuechstay;
GA_cue_stay.avg = (GA_erf.cuechstay.avg + GA_erf.cueenstay.avg) / 2;
GA_cue_switch = GA_erf.cuechswitch;
GA_cue_switch.avg = (GA_erf.cuechswitch.avg + GA_erf.cueenswitch.avg) / 2;

% then, calc the diff btwn the 2 conds
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_cue_stayvsswitch = ft_math(cfg, GA_cue_stay, GA_cue_switch);


% define parameters for plotting
start_time = stat.cfg.latency(1); % get the time window specified earlier in stat analysis
end_time = stat.cfg.latency(end);
timestep = 0.05; %(end_time - start_time) / 15; % length of time interval you want in each subplot (in seconds); 
                                                % alt: specify how many subplots you want (e.g. 15)
sampling_rate = 200; % we downsampled to 200Hz
sample_count = length(stat.time); % number of samples in MEG data (in the ERF time window)
j = [start_time : timestep : end_time];   % define the time interval (in seconds) for each subplot
m = [1 : timestep*sampling_rate : sample_count];  % corresponding sample indices in MEG data

% ensure stat.cfg.alpha (the alpha level we specified earlier in ft_timelockstatistics) still exists
if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.05; end; % if not, set it (new version corrects for 2-tailedness, so no need to use 0.025)

%{
if (length(stat.posclusters) == 0) % if no clusters were found at all, code below will throw error
    % so create a fake one (just to allow code below to run w/o error)
    stat.posclusters(1).prob = 1; 
    stat.posclusters(1).clusterstat = -9; 
    stat.posclusters(1).stddev = 0; 
    stat.posclusters(1).cirange = 0;
end
if (length(stat.negclusters) == 0) % do the same for neg clusters
    stat.negclusters(1).prob = 1; 
    stat.negclusters(1).clusterstat = -9; 
    stat.negclusters(1).stddev = 0; 
    stat.negclusters(1).cirange = 0;
end
%}

% get all p-values associated with the clusters
pos_cluster_pvals = [stat.posclusters(:).prob];
neg_cluster_pvals = [stat.negclusters(:).prob];
% find which clusters are significant, outputting their indices as held in stat.posclusters
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust); % I think stat.mask is simply combining pos & neg 
neg = ismember(stat.negclusterslabelmat, neg_signif_clust); % (i.e. stat.mask == pos | neg)

% Ensure the channels have the same order in the grand average and in the statistical output
% This might not be the case, because ft_math might shuffle the order  
[i1,i2] = match_str(GA_cue_stayvsswitch.label, stat.label);
% i1 holds a list of channel numbers in the grand averages
% i2 holds a list of channel numbers in the stat output

figure;  
for k = 1:length(j)-1; % create one subplot for each time interval
     subplot(3,5,k); % 3 * 5 = 15 subplots 
     
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   
     %cfg.zlim = [-5e-14 5e-14];  % set scaling (usually using automatic is ok) 
     pos_int = zeros(numel(GA_cue_stayvsswitch.label),1); % initialise the arrays with 0s
     neg_int = zeros(numel(GA_cue_stayvsswitch.label),1);
     pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2); % if a channel maintains significance thruout this time interval, then
     neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2); % we set this channel to 1 (to be highlighted)
     % not sure why it has to "maintain significance"; here I try with only requiring sig for half of time pts in this interval
     a = neg(i2, m(k):m(k+1));
     neg_int(i1) = sum(a, 2) > size(a, 2) / 2;
     
     sig_channels = find(pos_int | neg_int); % get indices of all significant channels
     if length(sig_channels) ~= 0 % if any sig channels found, report which channels these are
         fprintf(['In time interval [' num2str(cfg.xlim) '], these channels were significant:\n']);
         stat.label(sig_channels)
     end
     cfg.highlight = 'on';
     cfg.highlightchannel = sig_channels; % highlight these channels on topoplot
     cfg.highlightcolor = [255/255 192/255 203/255]; % pink colour

     cfg.comment = ['time = [' num2str(cfg.xlim) ']   ' strjoin(stat.label(sig_channels))]; % display time interval & names of sig channels
     %cfg.comment = 'auto'; % display date, xlim (time interval), zlim (amplitude range)
     cfg.commentpos = 'title';   
     %cfg.colorbar = 'yes'; % shows the scaling
     cfg.layout = lay;
     ft_topoplotER(cfg, GA_cue_stayvsswitch);
end  
%}


%% To plot the actual effect (i.e. average ERF of sig channels)

% SELECT which stat to plot & SPECIFY the relevant conds accordingly
stat = MixCost_interaction;
conds_to_plot = [1 3 4 6 7 9];
%stat = Main_Context;
%conds_to_plot = 1:9;%[1 2 4 5 7 8];

% load the relevant GA (avoid reloading if already exists)
if ~exist('GA_erf', 'var')
    load([ResultsFolder_thisrun 'GA_avg.mat']);
end

% read out the sig channels & time points
[row,col] = find(stat.mask);
sig_chans = unique(row)';
sig_samples = unique(col)';
start_time = stat.time(sig_samples(1));
end_time = stat.time(sig_samples(end));

% plot the average ERF over all sig channels
figure('Name', 'Avg ERF of sig channels');
cfg        = [];
cfg.title = ' '; % hide the display of channel names at the top
cfg.channel = stat.label(sig_chans);
cfg.baseline     = ERF_BASELINE; % makes no diff if we've already done baseline correction earlier
cfg.graphcolor   = cell2mat(colours(conds_to_plot)); 
cfg.linestyle    = lineTypes(conds_to_plot);
cfg.linewidth = 3;
cellarray = struct2cell(GA_erf);
ft_singleplotER(cfg, cellarray{conds_to_plot});
legend(eventnames_real(conds_to_plot), 'Location','northwest', 'FontSize',30);

xlim(PLOT_XLIM);
xlabel('Seconds');
ylabel('Tesla');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
% use alpha to set transparency 
alpha = 0.3;
patch(x,y,'black', 'FaceAlpha',alpha, 'HandleVisibility','off') % draw the shade 
    % (turn off HandleVisibility so it won't show up in the legends)
ylim(ylimits); % ensure ylim doesn't get expanded


% only plot the GFP if using axial ERF 
% (in the case of planar transformation, the entire timecourse is already +ve, so no need of GFP)
if ~PLANAR_TRANSFORM
    % if we want to plot the GFP of these channels, calculate that now
    cfg        = [];
    cfg.method = 'power';
    cfg.channel = stat.label(sig_chans);
    for j = 1:length(eventnames_real)
        GFP_Interaction.(eventnames_real{j}) = ft_globalmeanfield(cfg, GA_erf.(eventnames_real{j}));
    end
    
    % plot GFP
    figure('Name', 'GFP of sig channels'); hold on
    cfg       = [];
    cfg.title = ' '; % hide the display of channel names at the top
    cfg.graphcolor   = cell2mat(colours(conds_to_plot)); 
    cfg.linestyle    = lineTypes(conds_to_plot);
    cfg.linewidth = 3;
    ft_singleplotER(cfg, GFP_Interaction.NatStay, GFP_Interaction.NatSingle, GFP_Interaction.ArtStay, GFP_Interaction.ArtSingle, GFP_Interaction.BiStay, GFP_Interaction.BiSingle);
    legend(eventnames_real(conds_to_plot), 'Location','northwest', 'FontSize',30);
    
    xlim(PLOT_XLIM);
    xlabel('Seconds');
    ylabel('Tesla squared');
    set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
    box on; % draw a border around the figure

    % create shaded region indicating effect duration
    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
    y = [ylow ylow yhigh yhigh];
    % use alpha to set transparency 
    alpha = 0.3;
    patch(x,y,'black', 'FaceAlpha',alpha, 'HandleVisibility','off') % draw the shade 
        % (turn off HandleVisibility so it won't show up in the legends)
    ylim(ylimits); % ensure ylim doesn't get expanded
end


%% Ref code: copied directly from FT_compare_conditions.m
% see also: http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock#the_format_of_the_output
%{

neg_cluster_pvals=[];
pos_cluster_pvals=[];

figure;% % plot negative
cfg = [];
cfg.comment = 'no';
cfg.layout = layout;

cfg.xlim=latency;
cfg.highlight = 'off';

subplot(2,4,1)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,');']);
subplot(2,4,2)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,');']);
subplot(2,4,5)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_planar);']);
subplot(2,4,6)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,'_planar);']);

if isfield(eval([cond1,'_vs_',cond2,'_stat']),'posclusters') && ~isempty(eval([cond1,'_vs_',cond2,'_stat.posclusters']))
    eval(['pos_cluster_pvals = [',cond1,'_vs_',cond2,'_stat.posclusters(:).prob];']);
    pos_signif_clust = find(pos_cluster_pvals < 0.05);
    eval(['pos = ismember(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat, pos_signif_clust);']);
    eval(['poscluster_p=([',cond1,'_vs_',cond2,'_stat.posclusters.prob])']);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(pos);
    subplot(3,4,4)
    eval(['plot(GM_meg_',cond1,'.time,mean(GM_meg_',cond1,'.avg(logical(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat),:)))'])
    xlim([-0.2 0.5])
    xlabel('Time (s)')
    ylabel('Amplitude (fT)')
    hold on
    eval(['plot(GM_meg_',cond2,'.time,mean(GM_meg_',cond2,'.avg(logical(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat),:)))'])
else
    sprintf('no significant pos clusters')
end
subplot(2,4,3)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_vs_GM_meg_',cond2,');']);
cfg.highlight = 'off';
subplot(2,4,7)
eval(['ft_topoplotER(cfg, GM_meg_planar_',cond1,'_vs_GM_meg_planar_',cond2,');']);
set(gcf, 'Color', 'w');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 

figure;% % plot positive
cfg = [];
cfg.comment = 'no';
cfg.layout = layout;

cfg.xlim=latency;
cfg.highlight = 'off';

subplot(2,4,1)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,');']);
subplot(2,4,2)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,');']);
subplot(2,4,5)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_planar);']);
subplot(2,4,6)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,'_planar);']);

if isfield(eval([cond1,'_vs_',cond2,'_stat']),'negclusters') && ~isempty(eval([cond1,'_vs_',cond2,'_stat.negclusters']))
    eval(['neg_cluster_pvals = [',cond1,'_vs_',cond2,'_stat.negclusters(:).prob];']);
    neg_signif_clust = find(neg_cluster_pvals < 0.05);
    eval(['neg = ismember(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat, neg_signif_clust);']);
    eval(['negcluster_p=([',cond1,'_vs_',cond2,'_stat.negclusters.prob])']);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(neg);
    subplot(3,4,4)
    eval(['plot(GM_meg_',cond1,'.time,mean(GM_meg_',cond1,'.avg(logical(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat),:)))'])
    xlim([-0.2 0.5])
    xlabel('Time (s)')
    ylabel('Amplitude (fT)')
    hold on
    eval(['plot(GM_meg_',cond2,'.time,mean(GM_meg_',cond2,'.avg(logical(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat),:)))'])
else
    sprintf('no significant neg clusters')
end
subplot(2,4,3)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_vs_GM_meg_',cond2,');']);
cfg.highlight = 'off';
subplot(2,4,7)
eval(['ft_topoplotER(cfg, GM_meg_planar_',cond1,'_vs_GM_meg_planar_',cond2,');']);
set(gcf, 'Color', 'w');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
end

%}
