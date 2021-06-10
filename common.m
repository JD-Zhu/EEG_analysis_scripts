% the common #define section for all scripts
%
% Warning: all of the global vars below are CONSTANTS.
% UNDER NO CIRCUMSTANCES should their values be assigned/modified in other scripts.
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function [] = common()

    % specify all paths as absolute paths, to avoid any issues when we 'cd' into diff folders    
    global DataFolder; global ResultsFolder; global ResultsFolder_ROI; global ResultsFolder_Source;
    
    DataFolder = 'Z:\Analysis\Preprocess\SpinalCordInjury_CBD\EEG\DATA\'; % this directory should contain all the SubjectFolders
    ResultsFolder = 'Z:\Analysis\Preprocess\SpinalCordInjury_CBD\EEG\results_freq\'; % all subjects' freq analysis results will be stored here
    %DataFolder = [pwd '\\..\\DATA\']; % this directory should contain all the SubjectFolders
    %ResultsFolder = [pwd '\\..\\results_FREQ\\']; % all subjects' ERF results will be stored here
    %ResultsFolder_ROI = [pwd '\\..\\results_ROI\\']; % all subjects' ROI source-reconstruction results will be stored here
    %ResultsFolder_Source = [pwd '\\..\\results_SOURCE\\']; % all subjects' source localisation results will be stored here
    
    global filename_suffix; % select which pre-processing option: noPCA, reject components 1:3, or normal (reject components 1:5)
    % also need to change the last few lines in reject_response_component.m & load correct result files into the ResultsFolder
    filename_suffix = ''; % '_noPCA'; %'_rejectTop3'; %'';
    
    % actual EEG channels
    global EEG_chans;
    EEG_chans = [1:12 14:18 20:31 33:64]; % M1=13, M2=19, EOG=32
        
    % =================================================================
    
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
    
    % =================================================================
    
    % addpath to access custom functions in all subfolders
    addpath(genpath(pwd));
    
    % =================================================================
    
    % = Plot settings = %

    % time-course plots
    global PLOT_XLIM; global ERF_BASELINE; global ROI_BASELINE;
    PLOT_XLIM    = [-0.2 0.6];
    ERF_BASELINE = [-0.2 0];
    ROI_BASELINE = [-0.2 0];
    
    %global TFR_BASELINE;
    %TFR_BASELINE = [-0.75 0];
    
    % Plot shaded patch around time-course plots? 
    % Options: 'no', 'SEM', 'STDEV', 'CI_95'
    % (note: SEM < 95% CI < STDEV)
    global PLOT_SHADE;
    PLOT_SHADE = 'SEM';
    PLOT_SHADE = 'no'; % TEMP FIX - bounded_lines throws an error in stats_ROI!

    % Do we want to use a combination of diff colours & line types to
    % distinguish btwn conds? if no, we'll only use diff colours
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


    % toolbox to plot shaded boundary around each timecourse (diff paths for diff computers)
    addpath(genpath('C:\Users\Judy\Documents\MATLAB\kakearney-boundedline-pkg-50f7e4b'));
    addpath(genpath('C:\Users\43606024\Documents\MATLAB\kakearney-boundedline-pkg-50f7e4b'));
    
    % toolbox to save figure exactly as it appears on screen
    addpath(genpath('C:\Users\Judy\Documents\MATLAB\altmany-export_fig-9676767'));
    addpath(genpath('C:\Users\43606024\Documents\MATLAB\altmany-export_fig-9676767'));
    
    % SPM toolbox - for 3x3 ANOVA
    addpath('C:\Users\43606024\Documents\MATLAB\spm12');

end