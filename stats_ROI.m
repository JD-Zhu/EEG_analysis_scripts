%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_ROI.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Grand average & statistical analysis on reconstructed ROI activities
%
%    Q1: do stats across subjects?
%    A: yes. Plot grand ave first, just to get an idea of what effect is there. Then run stats.
%
%    Q2. use ft_timelockstatistics (exactly the same way as doing stats on erf)?
%    A: yes.
%    (but ft_timelockstats takes in erf structures, here we don't have the .avg field, just naked data)
%    A: either hack it into that format, or use EEGlab (statcond, std_stat).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;

% run the #define section
global ResultsFolder_ROI; % all subjects' ROI data are stored here

global eventnames_real; global colours_and_lineTypes; 
global colours; global lineTypes;
global PLOT_XLIM; global ROI_BASELINE;
global PLOT_SHADE; % for plotting shaded boundary on each time course
common();


% SELECT which set of single-subject ERFs to use
run_name = 'TSPCA10000_3'; % this should be a folder name inside the "Results_ROI" folder

% SELECT which beamformer results to use (fixed or free dipole orientation)
%run_suffix = '';         % for fixed orientation
run_suffix = '_freeori';  % for free orientation
%run_suffix = '_freeori_median-3';  % for free orientation, with smoothing applied on all single-subject ROI timecourses
    
ResultsFolder_ROI_thisrun = [ResultsFolder_ROI run_name run_suffix '\\'];


%% Read data

% find all .mat files in ResultsFolder_ROI_thisrun
files = dir([ResultsFolder_ROI_thisrun '*_ROI.mat']);

% each cycle reads in one '.mat' file (ie. one subject's ROI results)
for i = 1:length(files)
    filename = [ResultsFolder_ROI_thisrun files(i).name];
    load(filename);
    allSubjects_ROIs_bySubjects(i) = ROI_activity;
end

% get a list of all the ROI labels
ROIs_label = fieldnames(allSubjects_ROIs_bySubjects(1));

% reformat allSubjects_ROIs: Subject|ROI|condition -> ROI|condition|Subjects
for k = 1:length(ROIs_label)
    ROI_name = ROIs_label{k};
    allSubjects_ROIs.(ROI_name) = allSubjects_reformat(allSubjects_ROIs_bySubjects, ROI_name, eventnames_real);
%{
    for j = 1:length(eventnames_8) % 4 conditions in cue & 4 conditions in target (total 8)
       % if exist('allSubjects_ROI') % if var already exists, append to it
        %if ~isempty(['allSubjects_ROI.' ROIs_label{k} '.' (eventnames_8{j})]) % if var already exists, append to it
        if isfield(allSubjects_ROI, ROIs_label{k}) % if var already exists, append to it
            allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) = [allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) ROI_activity.(ROIs_label{k}).(eventnames_8{j})];
        else % if first time, simply assign to the first item
            allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) = ROI_activity.(ROIs_label{k}).(eventnames_8{j});
        end
    end
%}
end


%% Grand average across subjects

fprintf('\n= COMPUTING & PLOTTING CROSS-SUBJECT AVERAGES =\n');

% each cycle processes one ROI
for k = 1:length(ROIs_label)
    ROI_name = ROIs_label{k};
    fprintf(['\nROI: ' ROI_name '\n']);

    % compute grand average (across subjects) for each condition
    cfg = [];
    cfg.channel   = {'all'}; % there is only one channel (i.e. the virtual sensor for this ROI)
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    for j = 1:length(eventnames_real)
        cfg.keepindividual = 'no'; % average across subjects
        GA_avg.(eventnames_real{j}) = ft_timelockgrandaverage(cfg, allSubjects_ROIs.(ROI_name).(eventnames_real{j}){:}); 
        
        cfg.keepindividual = 'yes'; % do not average across subjects, keep the data for each individual subject
        GA_keepindi.(eventnames_real{j}) = ft_timelockgrandaverage(cfg, allSubjects_ROIs.(ROI_name).(eventnames_real{j}){:});

        % "{:}" means to use data from all elements of the variable
    end

    GA.(ROI_name) = GA_avg; % store it in the correct field
    GA_indi.(ROI_name) = GA_keepindi; % store it in the correct field
    
    % Plot the GAs
    %{
    figure('Name', ['GA in ' ROI_name]); hold on
        for j = 1:length(eventnames_real)
            plot(GA.(ROI_name).(eventnames_real{j}).time, GA.(ROI_name).(eventnames_real{j}).avg);
        end
    xlim(PLOT_XLIM);
    legend(eventnames_real);
    %}
end

% save the GA files
GA_output_file = [ResultsFolder_ROI_thisrun 'GA_avg.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA');
end
GA_output_file = [ResultsFolder_ROI_thisrun 'GA_individuals.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_indi');
end


%% Statistical analysis: prep for running ANOVA in SPM12
% Only need to run this section once & save the images
%{
fprintf('\n= STATS: ANOVA in SPM =\n');

SPM_temp_dir = [ResultsFolder_ROI_thisrun '\\SPM_images\\'];
mkdir(SPM_temp_dir);

% each cycle processes one ROI
for k = 1:length(ROIs_label)
    
    ROI_name = ROIs_label{k};
    fprintf(['\nROI: ' ROI_name '\n']);

    data = allSubjects_ROIs.(ROI_name); % data for the current ROI
    
    % Convert FT data to SPM format
    spm('defaults', 'eeg'); % make sure the path is correct
    
    conds_names = fieldnames(data);
    for j = 1:length(conds_names) % each cycle handles one cond
        for i = 1:length(data.(conds_names{j})) % each cycle handles one subject
            SPM_filename = [SPM_temp_dir ROI_name '_' conds_names{j} '_subj' num2str(i)];
            
            spm_eeg_ft2spm(data.(conds_names{j}){i}, SPM_filename);

            % change the type (or else it will throw an error)
            D = spm_eeg_load(SPM_filename);
            D = type(D, 'evoked');
            save([SPM_filename '_evoked.mat'], 'D');
            
            % Then convert to images using spm_eeg_convert2images
            S.D = [SPM_filename '_evoked.mat'];
            S.mode = 'time';
            %S.conditions = ; % default == all  % you can put multiple conditions in the same SPM object
            S.timewin = [-0.1 0.75]; % default == [-Inf Inf]  % we'll keep the whole time window, i.e. -100 ~ 750
            [images, outroot] = spm_eeg_convert2images(S);
        end
    end
end

%%% How to proceed %%%
% Now all the images have been saved in "SPM_Temp" (in all the "..._evoked" folders).
% You can now open SPM gui, and build your 3x3 model.
% Altenatively: to use the saved batch scripts, run "SPM_batch.m".

%}

%% Statistical analysis (to identify time interval of each effect, i.e. temporal clusters)

fprintf('\n= STATS: CLUSTER-BASED PERMUTATION TESTS =\n');

ft_defaults % reset paths (just in case we ran SPM12 before & ruined the paths)

% each cycle processes one ROI
for k = 1:length(ROIs_label)
    
    ROI_name = ROIs_label{k};
    fprintf(['\nROI: ' ROI_name '\n']);

    data = allSubjects_ROIs.(ROI_name); % data for the current ROI
    
    % set some config for the statistical test
    cfg = [];
    cfg.channel = {'all'}; % there is only one channel (i.e. the virtual sensor for this ROI)
    cfg.avgoverchan = 'yes'; % this is necessary (or else FT will ask for cfg.neighbours)
    
    cfg.latency = [-0.1 0.6]; % time interval over which the experimental 
                         % conditions must be compared (in seconds)
    cfg.avgovertime = 'no'; % if yes, this will average over the entire time window chosen in cfg.latency 
                            % (useful when you want to look at a particular component, e.g. to look at M100,
                            % cfg.latency = [0.08 0.12]; cfg.avgovertime = 'yes'; )

    %load([ResultsFolder_ROI_thisrun 'neighbours.mat']); % this is the sensor layout - it's the same for all subjects (even same across experiments). So just prepare once & save, then load here
    %cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

    cfg.method = 'montecarlo'; %'analytic';
    cfg.correctm = 'cluster'; %'no';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    %cfg.clusterstatistic = 'wcm'; cfg.wcm_weight = 1;    

    cfg.alpha = 0.05; % report all effects with p < 0.1
    cfg.numrandomization = 2000; % Rule of thumb: use 500, and double this number if it turns out 
        % that the p-value differs from the critical alpha-level (0.05 or 0.01) by less than 0.02

    numSubjects = length(files);
    within_design_1x2 = zeros(2, 2*numSubjects);
    within_design_1x2(1, :) = repmat(1:numSubjects, 1, 2);
    within_design_1x2(2, 1:numSubjects) = 1;
    within_design_1x2(2, numSubjects+1:2*numSubjects) = 2;

    within_design_1x3 = zeros(2, 3*numSubjects);
    within_design_1x3(1, :) = repmat(1:numSubjects, 1, 3);
    within_design_1x3(2, 1:numSubjects) = 1;
    within_design_1x3(2, numSubjects+1:2*numSubjects) = 2;
    within_design_1x3(2, 2*numSubjects+1:3*numSubjects) = 3;
    
    cfg.uvar  = 1; % row of design matrix that contains unit variable (in this case: subjects)
    cfg.ivar  = 2; % row of design matrix that contains independent variable (i.e. the conditions)

    
    %----- Run the statistical tests -----%
     
    % Use 'F test' for interaction & main effects (coz there are 3 levels in "context")
    cfg.statistic = 'ft_statfun_depsamplesFunivariate';
    cfg.design = within_design_1x3;
    cfg.tail = 1; % -1 = left, 1 = right
    cfg.clustertail = 1; % for F test, can only select right-sided tail
                         % https://github.com/fieldtrip/fieldtrip/blob/master/statfun/ft_statfun_depsamplesFunivariate.m

    % INTERACTION (i.e. calc sw$ in each context, then submit the 3 sw$ to F-test)
    fprintf('\n= Sw$ interaction (i.e. compare sw$ in 3 contexts using an F test) =\n');
    [timelock_SwCost_Nat, timelock_SwCost_Bi] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatStay, data.NatSwitch, data.BiStay, data.BiSwitch); %'2-1 vs 4-3'
    [timelock_SwCost_Art, ~] = combine_conds_for_T_test('fieldtrip', 'interaction', data.ArtStay, data.ArtSwitch, data.BiStay, data.BiSwitch);
    [SwCost_interaction.(ROI_name)] = ft_timelockstatistics(cfg, timelock_SwCost_Nat{:}, timelock_SwCost_Art{:}, timelock_SwCost_Bi{:});
    
    fprintf('\n= Mix$ interaction (i.e. compare mix$ in 3 contexts using an F test) =\n');
    [timelock_MixCost_Nat, timelock_MixCost_Bi] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatSingle, data.NatStay, data.BiSingle, data.BiStay);
    [timelock_MixCost_Art, ~] = combine_conds_for_T_test('fieldtrip', 'interaction', data.ArtSingle, data.ArtStay, data.BiSingle, data.BiStay);
    [MixCost_interaction.(ROI_name)] = ft_timelockstatistics(cfg, timelock_MixCost_Nat{:}, timelock_MixCost_Art{:}, timelock_MixCost_Bi{:});
       
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
    [Main_Context.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Nat{:}, timelock_Art{:}, timelock_Bi{:});

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
    [Main_Ttype.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Single{:}, timelock_Stay{:}, timelock_Switch{:});

    % Alternatively,
    % we can prob justify testing "switch effect" & "mixing effect" separately (no need to correct for MCP)
    fprintf('\nMain effect of switch: (t-test)\n');
    [Switch.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Stay{:}, timelock_Switch{:});
    fprintf('\nMain effect of mix: (t-test)\n');
    [Mix.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Single{:}, timelock_Stay{:});
    %fprintf('\nAlso testing single-vs-switch: (t-test)\n');
    %[Rubbish.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Single{:}, timelock_Switch{:});
    

    % PLANNED PAIRWISE COMPARISONS within each context 
    % (previously known as "SANITY CHECK" - did we find a switch cost in Bivalent context?)
    
    % Make sure we are using 2-tailed t-tests (see expla below in "UNPACKING main effects & interactions"):
    cfg.statistic = 'depsamplesT'; % t-test (i.e. for comparing 2 conds)
    cfg.design = within_design_1x2;
    cfg.tail = 0;
    cfg.clustertail = 0; % 2 tailed test
    cfg.correcttail = 'prob'; % correct for 2-tailedness
    
    fprintf('\n\n= Planned pairwise comparisons to assess sw$ & mix$ within each context\n');
    [stats_pairwise.Bi_sw.(ROI_name)] = ft_timelockstatistics(cfg, data.BiStay{:}, data.BiSwitch{:}); 
    [stats_pairwise.Bi_mix.(ROI_name)] = ft_timelockstatistics(cfg, data.BiSingle{:}, data.BiStay{:}); 
    [stats_pairwise.Nat_sw.(ROI_name)] = ft_timelockstatistics(cfg, data.NatStay{:}, data.NatSwitch{:}); 
    [stats_pairwise.Nat_mix.(ROI_name)] = ft_timelockstatistics(cfg, data.NatSingle{:}, data.NatStay{:}); 
    [stats_pairwise.Art_sw.(ROI_name)] = ft_timelockstatistics(cfg, data.ArtStay{:}, data.ArtSwitch{:}); 
    [stats_pairwise.Art_mix.(ROI_name)] = ft_timelockstatistics(cfg, data.ArtSingle{:}, data.ArtStay{:}); 
    
    contrasts = fieldnames(stats_pairwise);

    % write any sig effects to file
    fid = fopen([ResultsFolder_ROI_thisrun 'ROI_sanityCheck.txt'], 'at'); % open file for append
    for i = 1:length(contrasts)
        stat = stats_pairwise.(contrasts{i}).(ROI_name);
        if ~isempty(find(stat.mask))
            start_sample = find(stat.mask, 1, 'first');
            end_sample = find(stat.mask, 1, 'last');
            pvalue = stat.prob(start_sample); % p-value is the same for all time points in a cluster, so we just read it from the first time point
            start_time = stat.time(start_sample);
            end_time = stat.time(end_sample);
            fprintf(fid, '%s in %s, p = %.4f, between %.f~%.f ms (significant at samples %s).\n\n', ...
                contrasts{i}, ROI_name, pvalue, start_time*1000, end_time*1000, int2str(find(stat.mask))); % convert units to ms
        end
    end   
    fclose(fid);
   
end
    

%% UNPACKING main effects & interactions (only those found to be sig in F-tests above, specify the contrast & the ROI in which it was sig)
cfg.statistic = 'depsamplesT'; % t-test (i.e. for comparing 2 conds)
cfg.design = within_design_1x2;
%cfg.tail = -1; % -1 = left, 1 = right, 0 = 2-tailed
%cfg.clustertail = -1; % use left-sided tail, coz I always put the smaller cond first when calling ft_timelockstatistics()
% using 1-tailed t-test is generally frowned upon, so we change back to 2-tailed (and correct for it)
cfg.tail = 0; % -1 = left, 1 = right, 0 = 2-tailed
cfg.clustertail = 0; 
cfg.correcttail = 'prob'; % correct for 2-tailedness

% Unpack main effect of context
stat = Main_Context;
ROI_name = 'RACC';
cfg.latency = [0.125 0.150]; % duration of the cluster
cfg.avgovertime = 'yes';

data = allSubjects_ROIs.(ROI_name); % data for the current ROI
  
% run relevant section (during F-test above) to compute these:
% timelock_Nat, timelock_Art, timelock_Bi

[Context_nat_vs_bi.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Nat{:}, timelock_Bi{:});
[Context_art_vs_bi.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Art{:}, timelock_Bi{:});
[Context_nat_vs_art.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Nat{:}, timelock_Art{:});

length(find(Context_nat_vs_bi.(ROI_name).mask))
length(find(Context_art_vs_bi.(ROI_name).mask))
length(find(Context_nat_vs_art.(ROI_name).mask))


% We now use avgovertime/avgoverchan to unpack main effects & interactions,
% so the code below is obsolete.
%{ 
% To unpack the interactions, we compare the sw$ & mix$ for each pair of contexts (i.e. nat_vs_bi, art_vs_bi, nat_vs_art)
fprintf('\n\n= Sw$ Interaction - Unpacking (3 t-tests) =\n');
fprintf('\n  -> Nat vs Bi');
[SwCost_nat_vs_bi.(ROI_name)] = ft_timelockstatistics(cfg, timelock_SwCost_Nat{:}, timelock_SwCost_Bi{:});
fprintf('\n  -> Art vs Bi');
[SwCost_art_vs_bi.(ROI_name)] = ft_timelockstatistics(cfg, timelock_SwCost_Art{:}, timelock_SwCost_Bi{:});
fprintf('\n  -> Nat vs Art');
[SwCost_nat_vs_art.(ROI_name)] = ft_timelockstatistics(cfg, timelock_SwCost_Nat{:}, timelock_SwCost_Art{:});

fprintf('\n = Mix$ interaction - Unpacking (3 t-tests) =\n');
fprintf('\n  -> Nat vs Bi');
[MixCost_nat_vs_bi.(ROI_name)] = ft_timelockstatistics(cfg, timelock_MixCost_Nat{:}, timelock_MixCost_Bi{:});
fprintf('\n  -> Art vs Bi');
[MixCost_art_vs_bi.(ROI_name)] = ft_timelockstatistics(cfg, timelock_MixCost_Art{:}, timelock_MixCost_Bi{:});
fprintf('\n  -> Nat vs Art');
[MixCost_nat_vs_art.(ROI_name)] = ft_timelockstatistics(cfg, timelock_MixCost_Nat{:}, timelock_MixCost_Art{:});

% To unpack the main effect of context, we conduct 3 pairwise comparisons
fprintf('\n= Main effect of context - Unpacking (3 t-tests) =\n');
fprintf('\n  -> Nat vs Bi');
[Context_nat_vs_bi.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Nat{:}, timelock_Bi{:});
fprintf('\n  -> Art vs Bi');
[Context_art_vs_bi.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Art{:}, timelock_Bi{:});
fprintf('\n  -> Nat vs Art');
[Context_nat_vs_art.(ROI_name)] = ft_timelockstatistics(cfg, timelock_Nat{:}, timelock_Art{:});
%}

%save([ResultsFolder_ROI_thisrun 'stats_Interactions.mat'], 'SwCost_interaction', 'MixCost_interaction');
%save([ResultsFolder_ROI_thisrun 'stats_MainEffects.mat'], 'Main_Context', 'Main_Ttype', 'Switch', 'Mix');
%save([ResultsFolder_ROI_thisrun 'stats_pairwise.mat'], 'stats_pairwise');
%save([ResultsFolder_ROI_thisrun 'stats_Interactions_unpack.mat'], 'SwCost_nat_vs_bi', 'MixCost_nat_vs_bi', 'SwCost_art_vs_bi', 'MixCost_art_vs_bi', 'SwCost_nat_vs_art', 'MixCost_nat_vs_art');
%save([ResultsFolder_ROI_thisrun 'stats_MainEffects_unpack_avgovertime.mat'], 'Context_nat_vs_bi', 'Context_art_vs_bi', 'Context_nat_vs_art');


%% Find the effects & plot them
% Automatically check all the stats output & read out the time interval
% of each effect (from the stat.mask field)

close all;

stats = load([ResultsFolder_ROI_thisrun 'stats_Interactions.mat']); % select which stats output file to look at
load([ResultsFolder_ROI_thisrun 'GA_avg.mat']);
load([ResultsFolder_ROI_thisrun 'GA_individuals.mat']);

% make directory to store the output figures
mkdir([ResultsFolder_ROI_thisrun 'Figures\\non-sig\\']);

ROIs_names = fieldnames(GA); % get the list of ROI names

% Baseline correction b4 plotting. Notes: 
% 1. This is the right place to do baseline correction. We decided not to do it
% b4 running stats (i.e. at single-subject level) - see my email for expla
% 2. Only do baseline correction if we used "fixed dipole orientation" in beamformer.
% If we used "free dipole orientation" (i.e. entire timecourse is positive
% values), then we don't do baseline correction
if ~contains(ResultsFolder_ROI_thisrun, 'freeori')
    cfg = [];
    cfg.feedback = 'no';
    cfg.baseline = ROI_BASELINE; %[-0.1 0];
    for k = 1:length(ROIs_names) % each cycle handles one ROI
        ROI_name = ROIs_names{k};
        for j = 1:length(eventnames_real)
            GA.(ROI_name).(eventnames_real{j}) = ft_timelockbaseline(cfg, GA.(ROI_name).(eventnames_real{j})); 
        end
    end
end

% loop thru all the different contrasts & loop thru all ROIs in each,
% check whether the .mask field has any non-zero entries 
% (these are the effects & they are already cluster-corrected, so doesn't need to be consecutive 1s)
fprintf('\nThe following effects were detected:\n');
stats_names = fieldnames(stats);
for i = 1:length(stats_names) % each cycle handles one effect (e.g. cue_lang)
    stat_name = stats_names{i};
    ROIs_names = fieldnames(stats.(stat_name)); % get the list of ROI names
    
    for k = 1:length(ROIs_names) % each cycle handles one ROI
        ROI_name = ROIs_names{k};
                        
        % if the .mask contains any 1s, that's an effect (the rest are 0s)
        mask = stats.(stat_name).(ROI_name).mask; 
        effect = find(stats.(stat_name).(ROI_name).mask); 

        % do a GA plot for all contrasts that have an effect (will add gray box to show effect interval later)
        % also do a GA plot for all ROIs that don't have an effect (save in 'non-sig' folder); we don't need to plot all 3 contrasts for each ROI ('nat_vs_art', 'nat_vs_bi', 'art_vs_bi') - they are the same, just plot one
        if ( ~isempty(effect)) % || strcmp(stat_name(end-8:end), 'nat_vs_bi') ) % can't compare the whole word 'interaction' here, coz some stat_names (e.g. 'cue_lang') are a shorter string
            % GA plot
            figure('Name', [stat_name ' in ' ROI_name], 'Position', get(0, 'Screensize')); % make the figure full-screen size
            hold on;
            
            % grab the conds that are relevant for this stat_name
            type_of_cost = stat_name(1:2); % 'Sw' or 'Mix'
            %type_of_contrast = stat_name(end-9:end); % 'nat_vs_art' or 'nat_vs_bi'
            if strcmp(type_of_cost, 'Sw') % SwCost
                %{
                % if we only want to plot the 2 contexts that showed sig interaction:
                if strcmp(type_of_contrast, 'nat_vs_art')
                    conds_to_plot = [1 2 4 5];
                elseif strcmp(type_of_contrast, '_nat_vs_bi')
                    conds_to_plot = [1 2 7 8];
                else % 'art_vs_bi'
                    conds_to_plot = [4 5 7 8];
                end
                %}
                % if we want to plot all 3 contexts in same graph:
                conds_to_plot = [1 2 4 5 7 8];
            elseif strcmp(type_of_cost, 'Mi') % MixCost
                %{
                % if we only want to plot the 2 contexts that showed sig interaction:
                if strcmp(type_of_contrast, 'nat_vs_art')
                    conds_to_plot = [1 3 4 6];
                elseif strcmp(type_of_contrast, '_nat_vs_bi')
                    conds_to_plot = [1 3 7 9];
                else % 'art_vs_bi'
                    conds_to_plot = [4 6 7 9];
                end 
                %}
                % if we want to plot all 3 contexts in same graph:
                conds_to_plot = [1 3 4 6 7 9];                
            else % in all other cases (e.g. main effects, pairwise comparisons), plot all conds
                conds_to_plot = 1:9;
            end
            eventnames_subset = eventnames_real(conds_to_plot); 
            colours_subset = colours(conds_to_plot);
            lineTypes_subset = lineTypes(conds_to_plot);
                
            % each cycle plots 1 line (ie. 1 condition)
            for j = 1:length(eventnames_subset)
                if strcmp(PLOT_SHADE, 'no') % do not plot shaded boundary, just plot a single line                  
                    plot(GA.(ROI_name).(eventnames_subset{j}).time, GA.(ROI_name).(eventnames_subset{j}).avg, 'Color',colours_subset{j}, 'LineStyle',lineTypes_subset{j});
                else % calc the margin for shaded boundary (stdev / sem / CI) at every time point
                    allsubjects = GA_indi.(ROI_name).(eventnames_subset{j}).individual;
                    margin = calc_margin(allsubjects, PLOT_SHADE);

                    % plot time course with shaded boundary
                    boundedline(GA.(ROI_name).(eventnames_subset{j}).time, GA.(ROI_name).(eventnames_subset{j}).avg, margin(:), 'alpha', 'transparency',0.15, colours(j));                        
                end
            end

            xlim(PLOT_XLIM); 
            
            % set properties for axes, lines, and text
            xlabel('Seconds');
            ylabel('Ampere per square metre');
            set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
            box on; % draw a border around the figure

            % specify the legend manually (otherwise it will include
            % each shaded patch as an item too). For some reason,
            % the order of the lines are reversed when you grab them
            lines = findall(gcf, 'Type','line');
            legend(flip(lines(1:end)), ... % flip the order back to normal
              eventnames_subset, 'Location','northwest', 'FontSize',30);
            set(lines, 'Linewidth',3); % line thickness
                
            % reference code:
            %{
            p = findobj(gcf); % get the handles associated with the current figure
            
            allaxes = findall(p,'Type','axes');
            alllines = findall(p,'Type','line');
            alltext = findall(p,'Type','text');
            
            set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,'FontSize',14);
            set(alllines,'Linewidth',3);
            set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14);
            %}


            % if this GA plot is for a non-sig ROI, then save the figure into the 'non-sig' folder
            if isempty(effect)
                
                % for these plots, we want to delete the legends
                hl = findobj(gcf, 'type','legend');
                delete(hl);
                
                % maximise the figure before saving
                %set(gcf, 'Position', get(0, 'Screensize'));

                filename = [ROI_name '_' stat_name '.png'];
                %saveas(gcf, [ResultsFolder_ROI_thisrun 'Figures\\non-sig\\' filename]); % this fn does not maintain the aspect ratio, font size, etc
                export_fig(gcf, [ResultsFolder_ROI_thisrun 'Figures\\non-sig\\' filename]); % use this tool to save the figure exactly as shown on screen

            else % if there is any effect present, find all the clusters so we can output to console & mark on the plot

                % find the start & end of each cluster
                start_points = find(diff(mask) == 1); % transitions from 0 to 1
                end_points = find(diff(mask) == -1); % transitions from 1 to 0

                % check if we are missing the very first start point
                if (mask(1) == 1)
                    start_points = [1 start_points];
                end
                % check if we are missing the very last end point
                if (mask(end) == 1)
                    end_points = [end_points length(mask)];
                end
                % sanity check
                assert (length(start_points) == length(end_points));

                % increase the index for all start points by 1 (coz the "transition" found by diff 
                % is the position of the last '0' before it turns into '1')
                start_points = start_points + 1;


                % produce console output for each cluster & mark it on the GA plot
                fprintf('%s has an effect in %s:\n', ROI_name, stat_name);

                for cluster = 1:length(start_points) % start_points contains a list of starting points (one for each cluster)
                    start_sample = start_points(cluster);
                    end_sample = end_points(cluster);

                    % read out the necessary info
                    pvalue = stats.(stat_name).(ROI_name).prob(start_sample); % p-value is the same for all time points in a cluster, so we just read it from the first time point
                    start_time = stats.(stat_name).(ROI_name).time(start_sample);
                    end_time = stats.(stat_name).(ROI_name).time(end_sample); 

                    fprintf('    p = %.4f, between %.f~%.f ms (significant at samples %d to %d).\n', ...
                        pvalue, start_time*1000, end_time*1000, start_sample, end_sample); % convert units to ms

                    % mark the cluster interval on the GA plot
                    %line([start_time start_time], ylim, 'Color','black'); % plot a vertical line at start_time
                    %line([end_time end_time], ylim, 'Color','black'); % plot a vertical line at end_time

                    % create shaded region indicating effect duration
                    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
                    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
                    y = [ylow ylow yhigh yhigh];
                    % use alpha to set transparency 
                    if (pvalue < 0.05) % significant effect: use darker colour
                        alpha = 0.3;
                    else % marginal effect: use lighter colour
                        alpha = 0.05;
                    end
                    patch(x,y,'black', 'FaceAlpha',alpha, 'HandleVisibility','off') % draw the shade 
                        % (turn off HandleVisibility so it won't show up in the legends)
                    ylim(ylimits); % ensure ylim doesn't get expanded
                end

                % maximise the figure before saving
                %set(gcf, 'Position', get(0, 'Screensize'));
                
                % save the figure
                filename = [ROI_name '_' stat_name '.png'];
                %saveas(gcf, [ResultsFolder_ROI_thisrun 'Figures\\' filename]); % this fn does not maintain the aspect ratio, font size, etc
                export_fig(gcf, [ResultsFolder_ROI_thisrun 'Figures\\' filename]); % use this tool to save the figure exactly as shown on screen
                
                % old code to check multiple clusters
                %{
                % if the .mask contains any non-zero entries, that's an effect
                effect = find(stats.(stat_name).(ROI_name).mask); 
                if ~isempty(effect) % if there is an effect, we print it out & plot the ROI timecourse for each cond

                % check whether there are any jumps in the indices (i.e. multiple clusters)
                jump_positions = find(diff(effect) ~= 1);

                % find the start point for the first temporal cluster
                start_sample = effect(1);

                % find the end point for each cluster & print out this cluster
                for cluster = 1:length(jump_positions)

                    end_sample = effect(jump_positions(cluster));

                    % read out the necessary info
                    pvalue = stats.(stat_name).(ROI_name).prob(start_sample); % p-value is the same for all time points in a cluster, so we just read it from the first time point
                    start_time = stats.(stat_name).(ROI_name).time(start_sample);
                    end_time = stats.(stat_name).(ROI_name).time(end_sample); 

                    fprintf('%s has an effect in %s (p = %.4f), between %.f~%.f ms (significant at samples %d to %d).\n', ROI_name, stat_name, pvalue, start_time*1000, end_time*1000, start_sample, end_sample); % convert units to ms

                    % start point for the next cluster
                    start_sample = effect(jump_positions(cluster) + 1);
                end

                % find the end point for the final cluster & print out this cluster
                end_sample = effect(end); 
                % read out the necessary info
                pvalue = stats.(stat_name).(ROI_name).prob(start_sample); % p-value is the same for all time points in a cluster, so we just read it from the first time point
                start_time = stats.(stat_name).(ROI_name).time(start_sample);
                end_time = stats.(stat_name).(ROI_name).time(end_sample); 

                fprintf('%s has an effect in %s (p = %.4f), between %.f~%.f ms (significant at samples %d to %d).\n', ROI_name, stat_name, pvalue, start_time*1000, end_time*1000, start_sample, end_sample); % convert units to ms

                %}
            end
            
            hold off;
        end
    end
end
