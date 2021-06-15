% EEG frequency analysis
%
% Author: Judy Zhu (github.com/JD-Zhu)
%

%%
%close all
%clear all % disable this line if u want breakpoints to work

% run the #define section
global DataFolder; global ResultsFolder; 
% global EEG_chans; global colours;
common();


% find all subject folders containing raw EEG recording
SubjectIDs = dir([DataFolder '*_S1']);
%SubjectIDs = [dir([DataFolder 'A*']); dir([DataFolder 'B*'])];
%SubjectIDs = SubjectIDs([2 3 6]); % only process selected subjects
%SubjectIDs([2 13 25]) = []; % remove certain subjects from the list

SubjectIDs = {SubjectIDs.name}; % extract the names into a cell array
SubjectIDs = {'9007_S1'}; % or manually specify which subjects to process


% === Settings ===

% Please adjust as required:

% > create a name for this run (this will create a separate output & Figures folder)
%run_name = 'noICA';
run_name = 'offlineHPF';
output_name = ['output\\' run_name '\\']; % location to save intermediate output files inside each SubjectFolder
ResultsFolder_thisrun = [ResultsFolder run_name '\\']; % results for all subjects
mkdir(ResultsFolder_thisrun);

% > name of confile in the SubjectFolder (can use wildcards)
confile_name = '*.eeg';

% > which steps to run?
DO_HPF = true;
DO_ICA = true; % if we tested human subjects (i.e. not "dry run"), set this to true
RUN_ICA_ON_1HZ_FILTERED_DATA = true; % for MEG (a lot more channels), we prob don't need to apply 1Hz HPF before running ICA
                                     % for EEG, this step is recommended, otherwise ICA will just detect all the slow drifts & nothing useful
                                     % (https://www.youtube.com/watch?v=2hrYEYSycGI    https://jinjeon.me/post/eeg-advanced/)
                                     % I tried it - ICA decomposition was indeed poor quality if we don't apply 1Hz HPF first
DO_BEH_CHECK = false; % if subjects produced beh responses, set this to true
DO_PCA = false; % if subjects produced vocal responses, set this to true

% when running many subjects in one batch, process all auto steps until the first manual step
RUN_UP_TO_BEFORE_MANUAL_ARTEFACT = false;   % auto processing before 1st manual step
RUN_UP_TO_AFTER_MANUAL_ARTEFACT = false;    % perform 1st manual step (mark artefact)
RUN_UP_TO_ICA = false;                      % auto processing before 2nd manual step (ICA component analysis)
RUN_UP_TO_ICA_REJECTION = false;             % perform 2nd manual step (select ICA comps to reject)

% > other options:
CHANNEL_REPAIR = true; % repair bad/rejected channels? set to true for EEG data, as channel rejection leads to unbalanced offline reref
%CALC_UNCLEANED_ERF = false; % calculate uncleaned erf? (for quality check of response-component rejection)


% =================

% set filenames for saving the output from each stage (so that we don't have to rerun the whole thing from beginning every time)
S1_output_filename = 'S1_preprocessed_data.mat'; % Stage 1 output (stored inside each Subject folder)
%S2_output_filename = 'S2_after_visual_rejection.mat'; % Stage 2 output (stored inside each Subject folder)
S3_output_filename = ['.mat']; % Final output (stored in ResultsFolder for all subjects)


% load layout & neighbours
%load('easycapM11.mat'); % easycap doesn't have PO5 & PO6
load('lay_AntNeuro64.mat'); % use our custom-made layout & neighbours
load('neighbours_AntNeuro64.mat');
load('all_labels_AntNeuro64.mat'); % list of 61 real channels (i.e. excluding M1 M2 EOG)
%figure; ft_plot_layout(lay);
        


%% Stage 1: preprocessing

for i = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    
    % create output folder for all the save files
    output_path = [SubjectFolder output_name];
    mkdir(output_path);


    S1_output_file = [output_path S1_output_filename];
    
    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S1_output_file, 'file') ~= 2)          
        
        %{
        Preprocessing steps:
        - filtering
        - manually mark bad sections & bad channels
        - remove M1 M2 EOG
        - ICA artefact rejection
        - channel repair (interpolation)
        - offline rereferencing
        % trigger-based trial definition (i.e. epoching)
        %}
        
        
        % get the eeg data file name 
        files = dir([SubjectFolder confile_name]);
        rawfile = [SubjectFolder files(1).name];

        
        % >>>
        % Step 1: filtering
        output_file = [output_path 'filtered.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)   
            [alldata] = filtering(rawfile, DO_HPF);
            
            % save now, because the 0.01Hz HPF TAKES FOREVER to run!
            if (DO_HPF) % to save disk space, only save if we did HPF
                save(output_file, 'alldata', '-v7.3');
            end
        else
            load(output_file);
        end

        % >>>
        % Step 2: detrend - no longer used
        
        % >>>
        % Step 3: manually mark artefact sections
        output_file = [output_path 'arft.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)   
            % Print out SubjectID so we know which subject we are working on
            fprintf(['\nCURRENT SUBJECT: ' SubjectID '\n\n']); 

            [arft] = mark_artefact(alldata);
            save(output_file, 'arft', '-v7.3');
        else
            load(output_file);
        end

        % reject the manually marked artefact
        arft.artfctdef.reject = 'partial'; 
                                      % 'nan'; 
                                      % This fills those sections with NaNs.
                                      % You will then need to write code to
                                      % remove any trials containing NaNs
                                      % before computing ERF.
        alldata = ft_rejectartifact(arft, alldata);

        
        % >>>
        % Step 4: manually mark bad channels & reject them
        output_file = [output_path 'selChLabel.mat'];
                
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)  
            
            % First remove these channels: M1, M2, EOG (these should not contain 
            % brain activity, so we don't want them in the avg reference or the ICA!)
            % Doing this here (before manually rejecting bad channels)
            % because the indices would be different after that
            cfg         = [];
            cfg.channel = all_labels; %EEG_chans; %[1:12 14:18 20:31 33:64]; % M1=13, M2=19, EOG=32
            alldata = ft_selectdata(cfg, alldata);       
 
            [selChLabel] = reject_bad_channels(alldata);
            save(output_file, 'selChLabel', '-v7.3');
        else
            load(output_file);
        end

        % remove the bad channels
        cfg                         = [];
        cfg.channel                 = selChLabel;
        alldata                     = ft_selectdata(cfg, alldata);
        
        % If running in batch, skip to next subject now
        if (RUN_UP_TO_AFTER_MANUAL_ARTEFACT)
            continue;
        end        

        
        % >>>
        % Step 5: ICA
        if (DO_ICA)
            output_file = [output_path 'afterICA.mat'];
        
            % if haven't already processed this before, do it now & save a copy
            if (exist(output_file, 'file') ~= 2)  

                % Run ICA to identify components (if haven't done this yet)
                output_file_ICA = [output_path 'ICA_comps.mat'];        
                if (exist(output_file_ICA, 'file') ~= 2)    
                    
                    if (RUN_ICA_ON_1HZ_FILTERED_DATA) % apply 1Hz HPF before running ICA
                        [comp] = ICA_run(true, rawfile, arft, selChLabel);
                    else % directly run ICA without applying 1Hz HPF
                        [comp] = ICA_run(false, alldata);
                    end

                    % Immediately save a copy, as ICA takes a long time to run
                    save(output_file_ICA, 'comp', '-v7.3');
                else
                    load(output_file_ICA);
                end

                % If running in batch, skip to next subject now
                if (RUN_UP_TO_ICA)
                    continue;
                end


                % Select which ICA components to reject
                % set filename for diary (to record the components selected)
                set(0,'DiaryFile', [output_path 'ICA_log.txt']);

                [alldata] = ICA_reject_comps(alldata, comp, lay);

                save(output_file, 'alldata', '-v7.3');
            else
                load(output_file);
            end      
        end
        
        % If running in batch, skip to next subject now
        if (RUN_UP_TO_ICA_REJECTION)
            continue;
        end
        

        % >>>
        % Step 6: offline rereferencing
        
        % Interpolate rejected channels here,
        % otherwise it would make the offline reref unbalanced. 
        % This also avoids having too few channels in the GA (only 
        % those channels present in all subjects will be kept in the GA)
        if (CHANNEL_REPAIR)
            % add "elec" field to the data struct (needed for channel repair)
            %elec = ft_read_sens(rawfile, 'senstype','eeg', 'fileformat','easycap_txt');
            load('elec.mat'); % just load the version we have already made
            alldata.elec = elec;
            
            alldata = repair_bad_channels(alldata, neighbours, all_labels);
        end

        % re-reference using avg of all channels
        % https://www.fieldtriptoolbox.org/example/rereference/
        % https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_cleaning/
        cfg = [];
        cfg.reref      = 'yes';
        cfg.implicitref = 'CPz'; % add the online ref channel back into the data
        cfg.refchannel = 'all'; % which channels to use for offline reref
        cfg.refmethod  = 'avg';
        alldata = ft_preprocessing(cfg, alldata);
        
        
        % >>>
        % Step 7: downsample the data for saving
        %all_blocks.time(1:end) = all_blocks.time(1); % this avoids numeric round off issues in the time axes upon resampling
        cfg            = [];
        cfg.resamplefs = 200; % sampling freq was 1000Hz, best to use a divisor of it (200Hz is commonly used)
        cfg.detrend    = 'no';
        all_blocks     = ft_resampledata(cfg, alldata);

        % SAVE preprocessed data - takes a while!!
        %save(S1_output_file, 'all_blocks', 'trialinfo_b', '-v7.3');
        save(S1_output_file, 'all_blocks', '-v7.3');
    end
end


%%  Stage 2: trial exclusions
%{
for k = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(k));
    SubjectFolder = [DataFolder SubjectID '\\'];

    output_path = [SubjectFolder output_name];
    S2_output_file = [output_path S2_output_filename];

    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S2_output_file, 'file') ~= 2)    

        load([output_path S1_output_filename]);
   
        
        % divide up the master event list, to create 1 list for each cond
        %events_allBlocks = identify_event_types(SubjectID, trialinfo_b);
        
        % in each list, remove the indices corresponding to error trials
        if (DO_BEH_CHECK)
            events_allBlocks = exclude_beh_errors(SubjectID, events_allBlocks);
        end
 
    
        % === PCA artefact removal ===
        % remove mouth-movement artefact by extracting main components from the "response" epochs
        % and projecting these out of all trials

        if (DO_PCA)
            [all_blocks_clean, response_comp] = remove_artefact_PCA(all_blocks, events_allBlocks, lay, 'response');
            %TODO: Based on B11-pilot, it seems that removing top 5 comps is too
            % much. The pairwise comparison plots show that this removes the N1 peak etc
            % maybe change to 1:3 comps or sth, be more conservative.
        
        else % if not doing ICA, just update the variable name
            all_blocks_clean = all_blocks;
        end

        
        % === Reject Outlier Trials ===

        % Print out SubjectID so we know which subject we are working on
        fprintf(['\nCURRENT SUBJECT: ' SubjectID '\n\n']); 

        % Display visual trial summary to reject outlier trials
        cfg              = [];
        cfg.feedback     = 'no'; % suppress console output (so that it's easy to find the SubjectID we printed out above)
        cfg.method       = 'summary';
        cfg.metric       = 'zvalue'; % default is 'var'
        cfg.keepchannel  = 'no';
        cfg.keeptrial    = 'nan'; % we keep the rejected trials as 'NaN' here,
            % because if we remove them, that will change the indices of all subsequent trials,
            % which will no longer match the indices we are using in events_allBlocks
        all_blocks_clean = ft_rejectvisual(cfg, all_blocks_clean);
        
    
        %save([output_path S2_output_filename], 'all_blocks_clean', 'events_allBlocks', '-v7.3'); 
        save([output_path S2_output_filename], 'all_blocks_clean', '-v7.3'); 
        % 'all_blocks' was not changed in Stage 2, so don't need to save again
    end
end
%}


%% Stage 3: freq analysis

for i = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    
    output_path = [SubjectFolder output_name];
    S3_output_file = [ResultsFolder_thisrun SubjectID S3_output_filename];

    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S3_output_file, 'file') ~= 2)

        %load([output_path S2_output_filename]);
        load([output_path S1_output_filename]);
        
        
        %Q: Should we cut here or not? cutting seems to promote the alpha 
        % band power a lot more (is this more normal or less normal?)
        % DavidM: unless your data is very long (i.e. too computationally
        % intensive), it's actually better to use the whole continuous data
        %A: We want to look at infraslow freq band, so can't cut it!
        
        % cut the data into 4-second segments 
        %{
        % There are diff options for cutting (length, overlap, tapers):
        % https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_freq/
        cfg = [];
        cfg.length  = 4;
        cfg.overlap = 0;
        all_blocks_segmented = ft_redefinetrial(cfg, all_blocks);
        %}

        cfg = [];
        cfg.output  = 'pow';
        cfg.channel = 'all';
        cfg.method  = 'mtmfft';
        cfg.taper   = 'boxcar';
        cfg.foi     = 0:0.005:30; % 1 / cfg1.length = 0.25  
                     % (the longer the segments, the more reso we can have here)
        freq         = ft_freqanalysis(cfg, all_blocks);

        % unsquaring (i.e. back to uV unit)
        %bu = freq.powspctrm;
        %freq.powspctrm = sqrt(freq.powspctrm);
        
        % plot power spectrum of a particular channel
        %{
        figure;
        plot(freq.freq, freq.powspctrm(12,:))
        xlabel('Frequency (Hz)');
        ylabel('absolute power (uV^2)');
        xlim([1 30]);
        %}

        % where to save the figures
        save_location = [SubjectFolder 'Figures\\' run_name '\\'];
        mkdir(save_location);
        
        % plot power spectrum for all channels (overlay)
        figure; hold on;
        for chan = 1:length(freq.label)
            plot(freq.freq, freq.powspctrm(chan,:))
        end
        xlim([1 30]);
        xlabel('Frequency (Hz)');
        ylabel('Absolute power (uV^2)');
        hold off;
        
        export_fig(gcf, [save_location 'powspctrm_allchans.png']); % use this tool to save the figure exactly as shown on screen

        % plot avg of all channels (log transformed)
        figure; plot(freq.freq, mean(log(freq.powspctrm)));
        xlim([0.01 30]);
        xlabel('Frequency (Hz)');
        ylabel('Power (log[uV^2]');
        
        export_fig(gcf, [save_location 'powspctrm_avg.png']); % use this tool to save the figure exactly as shown on screen
        
        %{
        figure; 
        cfg = [];
        cfg.layout = lay;
        ft_multiplotER(cfg, freq);
        %}

        % topoplot for each freq band
        plot_TFR_topo(freq, lay, 'infraslow', [0.03 0.06], save_location)
        plot_TFR_topo(freq, lay, 'theta', [4 8], save_location)
        plot_TFR_topo(freq, lay, 'alpha', [9 12], save_location)
        plot_TFR_topo(freq, lay, 'beta', [13 25], save_location)
        
        
        % SAVE all relevant variables from the workspace
        save(S3_output_file, 'SubjectFolder', 'run_name', 'freq', '-v7.3');       
    else
        load(S3_output_file);
    end
end
