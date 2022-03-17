% Settings for Episodic Migraine project
%
% this acts as the common #define section for all scripts
%
% Warning: all of the global vars below are treated as CONSTANTS.
% UNDER NO CIRCUMSTANCES should their values be assigned/modified in any other scripts.
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function [] = common()

    % PLEASE SPECIFY THE SETTINGS FOR YOUR PROJECT BELOW %
    
    % 1. Location for your EEG data analysis 
    % (use absolute paths, to avoid any issues when we 'cd' into diff folders)
    ProjectFolder = 'Z:\Analysis\Judy\EpisodicMigraine\';
    
    global SUBJ_GROUP;
    SUBJ_GROUP = 'controls'; % Options: 'patients', 'controls'    
    
    global DataFolder;   
    DataFolder = [ProjectFolder 'data\' SUBJ_GROUP '\']; % this directory should contain all the SubjectFolders
    ResultsFolder = [ProjectFolder 'results\' SUBJ_GROUP '\']; % all subjects' freq analysis results will be stored here
    ResultsFolder_conn = [ProjectFolder 'results_conn\' SUBJ_GROUP '\']; % all subjects' connectivity results will be stored here

    
    % 2. Specify a list of subjects to analyse
    global SubjectIDs;
    
    % Option 1: find all subject folders inside DataFolder
    SubjectIDs = dir([DataFolder 'Subject*']);
    % you can modify the list in a few ways:
    %SubjectIDs = [dir([DataFolder 'A*']); dir([DataFolder 'B*'])]; % combine two lists
    %SubjectIDs = SubjectIDs([2 3 6]); % only process selected subjects
    %SubjectIDs([2 13 25]) = []; % remove certain subjects from the list

    SubjectIDs = {SubjectIDs.name}; % extract the names into a cell array

    % Option 2: manually specify the subject codes
    migraineurs_12 = {'Subject_500', 'Subject_548', 'Subject_583', 'Subject_661', ...
            'Subject_664', 'Subject_671', 'Subject_673', 'Subject_677', ...
            'Subject_680', 'Subject_681', 'Subject_696', 'Subject_800'};
    migraineurs_new5 = {'Subject_205', 'Subject_207', 'Subject_208', 'Subject_209', 'Subject_210'};
    controls_12 = {'Subject_101', 'Subject_495', 'Subject_622', 'Subject_624', ...
            'Subject_634', 'Subject_642', 'Subject_675', 'Subject_690', ...
            'Subject_809', 'Subject_885', 'Subject_886', 'Subject_891'};
    controls_13_remaining = {'Subject_608', 'Subject_610', 'Subject_613', 'Subject_623', ...
                    'Subject_629', 'Subject_631', 'Subject_640', 'Subject_645', ...
                    'Subject_682', 'Subject_804', 'Subject_808', 'Subject_844', 'Subject_846'};
    controls_4new = {'Subject_251', 'Subject_252', 'Subject_253', 'Subject_254'};

    if strcmp(SUBJ_GROUP, 'patients')
        SubjectIDs = [migraineurs_12 migraineurs_new5];
    elseif strcmp(SUBJ_GROUP, 'controls')
        SubjectIDs = [controls_12 controls_13_remaining controls_4new];
    end

    % or process these new subjects only
    %SubjectIDs = {'Subject_608'}; %846, 640, 631, 629
    
    
    % 3. Which EEG system did you use to collect the data?
    EEG_system = 'NeuroPrax32'; % Options: 'AntNeuro64', 'NeuroPrax32'
    
    % these custom layout files for each system were made using prepare_layout_and_neighbours.m
    global LAYOUT_FILE; global NEIGHBOURS_FILE; global ALL_LABELS_FILE; global ELEC_FILE;
    LAYOUT_FILE = ['lay_' EEG_system '.mat'];
    NEIGHBOURS_FILE = ['neighbours_' EEG_system '.mat'];
    ALL_LABELS_FILE = ['all_labels_' EEG_system '.mat'];
    ELEC_FILE = ['elec_' EEG_system '.mat'];

    global CONFILE_NAME; % name of the raw data file in each subject folder (can use wildcards)
    global ONLINE_REF; % specify the online ref channel for this EEG system (if you want to add it back into the data during reref)
    if strcmp(EEG_system, 'AntNeuro64')
        CONFILE_NAME = '*.eeg';
        ONLINE_REF = 'CPz'; 
    elseif strcmp(EEG_system, 'NeuroPrax32')
        CONFILE_NAME = '*\*.edf';
        ONLINE_REF = ''; 
    end

    
    % 4. Options for preprocessing & single-subject analysis (main.m)
    
    % (1) create a name for this run (this will create a separate output & Figures folder)
    global run_name;
    run_name = '_EC_LPF30'; % '_EO';
    
    global file_suffix; % only used in main_migraine.m right now - can we get rid of this?
    file_suffix = ''; %'_minReject': only reject a noisy chan if it's utterly crazy - keep where possible (note: all flat channels must still be rejected)


    % (2) offline rereferencing using "average reference" or "linked mastoid"?
    global REREF; 
    REREF = 'AR'; % we don't have M1 M2 for these data, so just use avg ref

    if strcmp(REREF, 'LM')
       run_name = [run_name '_LMref']; 
    end
    
    
    % (3) for connectivity analysis ONLY - apply surface Laplacian to deal with volumn conduction issue?
    global APPLY_SL; 
    APPLY_SL = true;

    
    % (4) include infra-slow oscillations in the analysis?
    global ANALYSE_ISO;
    ANALYSE_ISO = false; % ISO range was already filtered out in these data
    
    
    % (5) which preprocessing steps to run?
    global DO_HPF; global FILTERS; global PLOT_CHANNEL_SPECTRA; 
    global DO_ICA; global FILTER_AGAIN_BEFORE_ICA; global FILTERS_for_ICA; 
    global CHANNEL_REPAIR; global DO_BEH_CHECK; global DO_PCA;
    global DOWNSAMPLE; % set this to 0 if no downsampling needed (e.g. for chronic migraine data acquired at 125Hz sampling rate)
    DO_HPF = true;
    FILTERS = [1 2 35 10]; % HPF 1+-1Hz; LPF 35+-5Hz
    PLOT_CHANNEL_SPECTRA = false; % during initial data inspection, plot channel spectra to help with determining bad channels?
                                  % Note: channel spectra is plotted on raw data (i.e. without filtering)
    DO_ICA = false; % short 5min recording, eyes-closed resting state
    FILTER_AGAIN_BEFORE_ICA = false; % this dataset is already HPF'd at 1Hz, no need to filter again before ICA
    %HPF_for_ICA = [1 2 FILTERS(3) FILTERS(4)]; % n/a
    CHANNEL_REPAIR = true; % interpolate rejected channels? 
                           % must set to true if using "average reference", as channel rejection leads to unbalanced reref
    DOWNSAMPLE = 250; % old episodic migraine data: sampling rate was 500Hz, best to use a divisor of it (~200Hz is commonly used)
    DO_BEH_CHECK = false; % if subjects produced beh responses, set this to true
    DO_PCA = false; % if subjects produced vocal responses, set this to true

    % when running many subjects in one batch, process all auto steps until the next manual step
    global RUN_UP_TO_BEFORE_MANUAL_ARTEFACT; global RUN_UP_TO_AFTER_MANUAL_ARTEFACT; 
    global RUN_UP_TO_ICA; global RUN_UP_TO_ICA_REJECTION; global BROWSING_WITHOUT_SAVE;
    RUN_UP_TO_BEFORE_MANUAL_ARTEFACT = false;   % auto processing before 1st manual step
    RUN_UP_TO_AFTER_MANUAL_ARTEFACT = false;    % perform 1st manual step (mark artefact & reject bad channels)
    RUN_UP_TO_ICA = false;                      % auto processing before 2nd manual step (ICA component analysis)
    RUN_UP_TO_ICA_REJECTION = false;            % perform 2nd manual step (select ICA comps to reject)
    BROWSING_WITHOUT_SAVE = false;              % browse filtered data - do not save arft & selChLabels

    % locations to save the output                   
    temp_name = run_name(2:end);
    if isempty(temp_name)
        temp_name = 'full';
    end

    global output_name; % for intermediate output files during preprocessing (so that we don't have to rerun the whole thing from beginning every time)
    global ResultsFolder_thisrun; global ResultsFolder_conn_thisrun; % results for all subjects
    output_name = ['output' run_name '\\']; % TODO (future): separate these from the DataFolder - put them in a separate "preprocessed" folder at the top level
    ResultsFolder_thisrun = [ResultsFolder temp_name '\\'];
    ResultsFolder_conn_thisrun = [ResultsFolder_conn temp_name '\\'];
    
    % filenames for saving the intermediate output from each stage of preprocessing
    global S1_output_filename; global S3_output_filename; global S4_output_filename;
    S1_output_filename = ['S1_preprocessed_data' file_suffix '.mat']; % Stage 1 output (stored inside each Subject folder under output_name)
    %S2_output_filename = ['S2_after_visual_rejection' file_suffix '.mat']; % Stage 2 output (stored inside each Subject folder under output_name)
    S3_output_filename = [file_suffix '.mat']; % Final output for freq analysis (stored in ResultsFolder for all subjects)
    S4_output_filename = [file_suffix '.mat']; % Final output for connectivity analysis (stored in ResultsFolder_conn for all subjects)

    
    % (6) special troubleshooting steps for old episodic migraine data only
    global EPISODIC_ONLY;
    EPISODIC_ONLY = true;
    
    
    % 5. Options for grand average & excel export (stats_FREQ.m)
    
    % (1) are we working with connectivity results here? 
    global is_conn; % this setting is not just for setting the correct folder below - it's being used in a number of places in stats_FREQ.m
    is_conn = true;
    
    
    % (2) specify a list of subjects to compute GA on
    
    % this can be subjects belonging to a particular sub-group (e.g. prodromes / postdromes / interictals),
    % or leave empty to use all subjects in the ResultsFolder (i.e. all patients / all controls)
    global SubjectIDs_GA;
    SubjectIDs_GA = [];
    % Episodic migraine proj - final set of 17 controls (age & gender matched to migraineurs)
    SubjectIDs_GA = {'Subject_101', 'Subject_251', 'Subject_252', 'Subject_253', 'Subject_254', 'Subject_495', 'Subject_610', 'Subject_622', 'Subject_623', 'Subject_634', 'Subject_642', 'Subject_675', 'Subject_690', 'Subject_809', 'Subject_844', 'Subject_885', 'Subject_891'};
    % Groups based on migraine phases:
    %SubjectIDs_GA = {'Subject_500', 'Subject_548', 'Subject_208'}; % prodrome
    %SubjectIDs_GA = {'Subject_583', 'Subject_673', 'Subject_680', 'Subject_205'}; % postdrome
    %SubjectIDs_GA = {'Subject_661', 'Subject_664', 'Subject_671', 'Subject_677', 'Subject_681', 'Subject_696', 'Subject_800', 'Subject_207', 'Subject_209', 'Subject_210'}; % interictal
    % Groups based on migraine frequency:
    %SubjectIDs_GA = {'Subject_677', 'Subject_681', 'Subject_696', 'Subject_800'}; % <1 day / month
    %SubjectIDs_GA = {'Subject_583', 'Subject_661', 'Subject_671'}; % 1-2 days / month
    %SubjectIDs_GA = {'Subject_500', 'Subject_548', 'Subject_664', 'Subject_673', 'Subject_680'}; % >3 days / month

    % if subject list is empty, then use all results files in the folder
    if isempty(SubjectIDs_GA)
        SubjectIDs_GA = dir([ResultsFolder_thisrun '*.mat']);
        SubjectIDs_GA = {SubjectIDs_GA.name}; % extract the names into a cell array
        SubjectIDs_GA = cellfun(@(x) x(1:end-4), SubjectIDs_GA, 'un', 0); % remove the '.mat' extension
    end

    
    % (3) exporting to excel
    global FREQ_FIELD; % for fixing up the "freq" field in results (for some reason the freqs are not whole numbers)
    global FREQS_TO_EXPORT; % if analysing ISO, only export certain freqs to excel (coz we computed 3001 freq points: 0:0.01:30)
    
    if ANALYSE_ISO
        FREQ_FIELD = 0:0.01:30; 
        FREQS_TO_EXPORT = [0.02:0.01:0.09 0.1:0.1:0.9 1:30];
    else
        FREQ_FIELD = 1:30;
        FREQS_TO_EXPORT = FREQ_FIELD;
    end

    
    % 6. Plot settings
    global PLOT_XLIM; %global ERF_BASELINE; global ROI_BASELINE;
    
    % for FREQ plots
    if ANALYSE_ISO
        PLOT_XLIM = [1 30]; % anything below 1Hz is way over powered (rendering the whole plot unviewable)
    else
        PLOT_XLIM = [2 30]; % for old episodic migraine data, we are interested in 2-30Hz (everything below 2Hz was already filtered out?)
    end
    
    % for ERP plots
    %PLOT_XLIM    = [-0.2 0.6];
    %ERF_BASELINE = [-0.2 0];
    %ROI_BASELINE = [-0.2 0];
    
    %global TFR_BASELINE;
    %TFR_BASELINE = [-0.75 0];
    
    % Plot shaded patch around time-course plots? 
    % Options: 'no', 'SEM', 'STDEV', 'CI_95'
    % (note: SEM < 95% CI < STDEV)
    global PLOT_SHADE;
    %PLOT_SHADE = 'SEM';
    PLOT_SHADE = 'no'; % TEMP FIX - bounded_lines throws an error in stats_ROI!

    % Do we want to use a combination of diff colours & line types to
    % distinguish btwn conds? if no, we'll only use diff colours
    %{
    global colours_and_lineTypes;
    colours_and_lineTypes = true;
    
    numConds = length(eventnames_real); % total number of conds to plot
    numCategories = 3; % ('Bi','Nat','Art') 
                       % All ttypes in each context will be same colour.
                       % Note: this var is only applicable when
                       % colours_and_lineTypes is set to 'true'.
    
    
    % colours for time course plots (one colour for each condition):
    % need to specify manually because we plot each cond separately, and simply 
    % using default colourmap makes all lines the same colour when calling boundedline()
    global colours; global lineTypes;
    
    if (colours_and_lineTypes)
        colours = {'g','g','g','b','b','b','r','r','r'};
        %colours = repelem(colour_list, numConds/numCategories, 1); % repeat each colour 3 times: [1 1 1 2 2 2 3 3 3]
        
        lineTypes = repmat({'-', '--', ':'}, [1 numConds/numCategories]); % {'-', '--', ':', '-', '--', ':', '-', '--', ':'}
    else
        % only 6 unique colours below, need more colours (we have 9 conds)!
        %colours = ['b', 'r', 'g', 'k', 'y', 'm', 'b', 'r', 'g', 'k', 'y', 'm'];
        colour_list = distinguishable_colors(numConds);

        % include the seq twice, once for the lines, once for the shaded boundaries
        %colours = ['b', 'r', 'y', 'm', 'b', 'r', 'y', 'm']; 
        colours = [colour_list; colour_list];
        
        % set line type to default (solid line)
        lineTypes = repmat({'-'}, [1 numConds]); % set to '-' for all conds
        lineTypes = [lineTypes lineTypes]; % repeat again for the shaded boundaries?
    end
    %}
    

    % =================================================================
    
    % addpath to access custom functions in all subfolders
    addpath(genpath(pwd));
    
    % toolbox to plot shaded boundary around each timecourse
    %addpath(genpath('C:\Users\Judy\Documents\MATLAB\kakearney-boundedline-pkg-50f7e4b'));
    
    % toolbox to save figure exactly as it appears on screen
    %addpath(genpath('C:\Users\Judy\Documents\MATLAB\altmany-export_fig-9676767'));

        
    % =================================================================
    
    %{
    % Using photodetector to adjust for screen delay
    % -> use each individual trigger time or the median screen delay?
    % Options: 'median', 'individual_trial'
    global ADJUST_SCREEN_DELAY;
    ADJUST_SCREEN_DELAY = 'median'; % 'individual' is more valid than 'median'
    
    
    % trigger events (DO NOT change the order of this list)
    global eventcodes; global eventnames; global eventnames_real;
    eventcodes = {...
    {'NatStayC'},{'17'};{'NatStayE'},{'18'};{'NatSwitchC'},{'19'};{'NatSwitchE'},{'20'};{'NatSingleC'},{'21'};{'NatSingleE'},{'22'}; ...
    {'ArtStayC'},{'23'};{'ArtStayE'},{'24'};{'ArtSwitchC'},{'25'};{'ArtSwitchE'},{'26'};{'ArtSingleC'},{'27'};{'ArtSingleE'},{'28'}; ...
    {'BiStayC'},{'37'};{'BiStayE'},{'38'};{'BiSwitchC'},{'39'};{'BiSwitchE'},{'40'};{'BiSingleC'},{'41'};{'BiSingleE'},{'42'}; ...
    {'response'},{'30'}};
    eventnames = eventcodes(:,1); % extract a list of all event names
    eventnames = [eventnames{:}]; % convert into strings
    eventnames_real = eventnames(1:18); % exclude the 'response' event, which is not a real cond
    % for AEF pilot, we don't have "bivalent" conds, so only use the first 12 conds
    %eventnames = eventnames(1:12); 
    
    % do we want to collapse across Eng & Chn?
    global collapse_across_langs;
    collapse_across_langs = true;
    
    if (collapse_across_langs)
        eventnames_real = eventnames_real([2,4,6,8,10,12,14,16,18]); % only take 1 from every pair of eventnames
        eventnames_real  = cellfun(@(x) x(1:end-1), eventnames_real, 'UniformOutput', false); % strip the last char (indicating which lang)
    end    
    %eventnames_real = eventnames_real([1 2 3 6 9]); % partial exp (Natuni & singleLang only)
    %}
    
    % from exp 1:
    %{
    eventcodes = {{'cuechstay'},{'17'};{'cuechswitch'},{'19'};{'cueenstay'},{'21'};{'cueenswitch'},{'23'}; ...
        {'targetchstay'},{'18'};{'targetchswitch'},{'20'};{'targetenstay'},{'22'};{'targetenswitch'},{'24'};{'response'},{'30'}};
    eventnames = eventcodes(:,1); % extract a list of all event names
    eventnames = [eventnames{:}]; % convert into strings
    
    % for ease of reference to the conditions in cue window & target window
    global conds_cue; global conds_target;
    conds_cue = 1:4;
    conds_target = 5:8;
    eventnames_8 = eventnames([conds_cue conds_target]); % 8 actual event types
    %}

end