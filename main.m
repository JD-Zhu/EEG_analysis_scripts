% EEG frequency analysis
%
% Author: Judy Zhu (github.com/JD-Zhu)
%

%%
%close all
%clear all % disable this line if u want breakpoints to work

% run the #define section
global DataFolder; global SUBJ_GROUP; global SubjectIDs; global CONFILE_NAME; 
global LAYOUT_FILE; global NEIGHBOURS_FILE; global ALL_LABELS_FILE; global ELEC_FILE;
global ONLINE_REF; global REREF; global APPLY_SL; global run_name; global file_suffix;
global output_name; global ResultsFolder_thisrun; global ResultsFolder_conn_thisrun;
global S1_output_filename; global S3_output_filename; global S4_output_filename;
global DO_HPF; global FILTERS; global PLOT_CHANNEL_SPECTRA; 
global DO_ICA; global FILTER_AGAIN_BEFORE_ICA; global FILTERS_for_ICA; 
global CHANNEL_REPAIR; global DOWNSAMPLE; %global DO_BEH_CHECK; global DO_PCA;
global RUN_UP_TO_BEFORE_MANUAL_ARTEFACT; global RUN_UP_TO_AFTER_MANUAL_ARTEFACT; 
global RUN_UP_TO_ICA; global RUN_UP_TO_ICA_REJECTION; global BROWSING_WITHOUT_SAVE;
global ANALYSE_ISO; global FREQ_BANDS; global EPISODIC_ONLY;
% global colours;
common();

% location to save the results for all subjects
mkdir(ResultsFolder_thisrun);
mkdir(ResultsFolder_conn_thisrun);

% load layout & neighbours
%load('easycapM11.mat'); % easycap doesn't have PO5 & PO6
load(LAYOUT_FILE); % use our custom-made layout & neighbours
%figure; ft_plot_layout(lay);
load(NEIGHBOURS_FILE);
load(ALL_LABELS_FILE); % list of real EEG channels (i.e. excluding EOG & ref chans): 61 chans for AntNeuro, 27 chans for NeuroPrax
        

%% Stage 1: preprocessing

for i = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    
    % create output folder for all the save files
    output_path = [SubjectFolder output_name];
    mkdir(output_path);
    
    disp(['   SubjectFolder: ' SubjectFolder]);
    disp(['   output_path    ' output_path]);


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
        files = dir(fullfile(SubjectFolder, CONFILE_NAME));
        rawfile = fullfile(files(1).folder, files(1).name);

            
        % >>>
        % Step 1: read in the raw data
        %hdr = ft_read_header(rawfile);%, 'dataformat','yokogawa_con'); % read header file
         
        % ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
        cfg                      = [];
        cfg.trialfun             = 'ft_trialfun_general';
        cfg.datafile             = rawfile;
        cfg.headerfile           = [rawfile(1:end-3) 'vhdr'];
        cfg.trialdef.triallength = Inf;
        cfg.trialdef.ntrials     = 1; % read in all data as a single segment, coz filtering should be done on continuous data
        cfg = ft_definetrial(cfg);

        % https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_cleaning/
        cfg.demean     = 'yes';
        cfg.detrend    = 'no';
        cfg.continuous = 'yes';
        alldata = ft_preprocessing(cfg);
    
    
        % >>>
        % Step 2: filtering
        output_file = [output_path 'filtered.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)   
            [alldata] = filtering(alldata, DO_HPF, FILTERS);
            
            % save now, because the 0.01Hz HPF TAKES FOREVER to run!
            if (DO_HPF) % to save disk space, only save if we did HPF
                save(output_file, 'alldata', '-v7.3');
            end
        else
            load(output_file);
        end

        
        % If running in batch, skip to next subject now
        if (RUN_UP_TO_BEFORE_MANUAL_ARTEFACT)
            continue;
        end
        
        
        % >>>
        % Step 3: manually mark artefact sections
        output_file = [output_path 'arft.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)   
            % Print out SubjectID so we know which subject we are working on
            fprintf(['\nCURRENT SUBJECT: ' SubjectID '\n\n']); 

            [arft] = mark_artefact(alldata, [-64 64]);
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
        
        
        % if using linked mastoid ref, keep a copy of M1 & M2 here
        if strcmp(REREF, 'LM')
            cfg = [];
            cfg.channel = {'M1', 'M2'};
            M1M2_data = ft_selectdata(cfg, alldata);
            
            save([output_path 'M1M2_data.mat'], 'M1M2_data', '-v7.3');
        end
        
        
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
                    
                    if (FILTER_AGAIN_BEFORE_ICA) % apply another HPF before running ICA
                        [comp] = ICA_run(true, FILTERS_for_ICA, rawfile, arft, selChLabel);
                    else % directly run ICA without applying another HPF
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

                [alldata] = ICA_reject_comps(alldata, comp, lay, output_path);
                close all;

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
        % Step 6: channel repair
        
        % Interpolate rejected channels here, otherwise it would make 
        % the offline reref unbalanced (if using avg ref - see below). 
        % This also avoids having too few channels in the GA (only 
        % those channels present in all subjects will be kept in the GA)
        if (CHANNEL_REPAIR)
            % add "elec" field to the data struct (needed for channel repair)
            %elec = ft_read_sens(rawfile, 'senstype','eeg', 'fileformat','easycap_txt');
            load(ELEC_FILE); % just load the version we have already made
            alldata.elec = elec;
            
            alldata = repair_bad_channels(alldata, neighbours, all_labels);
        end
         
        % Adjust the order of channels so that they are consistent across ptps
        % (this step is not needed for FT, but if you ever export the
        % data/results, e.g. to Excel, having a diff channel order for each
        % ptp would be very confusing/error-prone)
        % https://mailman.science.ru.nl/pipermail/fieldtrip/2017-March/037333.html
        % https://mailman.science.ru.nl/pipermail/fieldtrip/2015-March/034890.html
        
        idx = arrayfun( @(x)(find(strcmp(alldata.label, x))), all_labels);
        % adjust chan order in the label field
        alldata.label = alldata.label(idx);
        % adjust chan order in the data field (note that there may be 
        % multiple "trials", e.g. if we manually marked artefact earlier)
        for n = 1:length(alldata.trial)
            tmp = alldata.trial{n};
            tmp = tmp(idx,:);
            alldata.trial{n} = tmp;
        end

        
        % >>>
        % Step 7: offline rereferencing
        
        % https://www.fieldtriptoolbox.org/example/rereference/
        % https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_cleaning/
        
        if strcmp(REREF, 'AR') % re-reference using avg of all channels
            cfg = [];
            cfg.reref      = 'yes';
            cfg.refchannel = 'all'; % which channels to use for offline reref
            cfg.refmethod  = 'avg';
            if ~isempty(ONLINE_REF)
                cfg.implicitref = ONLINE_REF; % add the online ref channel back into the data (will be filled with 0)
            end
            alldata = ft_preprocessing(cfg, alldata);
        elseif strcmp(REREF, 'LM') % re-reference using linked mastoid (i.e. avg of M1 & M2)  
            % add M1 & M2 back in first
            load([output_path 'M1M2_data.mat']);
            alldata = ft_appenddata([], alldata, M1M2_data);
            
            cfg = [];
            cfg.reref      = 'yes';
            cfg.refchannel = {'M1', 'M2'}; % which channels to use for offline reref
            cfg.refmethod  = 'avg';
            if ~isempty(ONLINE_REF)
                cfg.implicitref = ONLINE_REF; % add the online ref channel back into the data (will be filled with 0)
            end
            alldata = ft_preprocessing(cfg, alldata);
            
            % remove M1 & M2
            cfg         = [];
            cfg.channel = {'all', '-M1', '-M2'};
            alldata = ft_selectdata(cfg, alldata);   
        end
        
        all_blocks = alldata;
        
        
        % >>>
        % Step 8: downsample the data for saving
        if DOWNSAMPLE ~= 0
            %all_blocks.time(1:end) = all_blocks.time(1); % this avoids numeric round off issues in the time axes upon resampling
            cfg            = [];
            cfg.resamplefs = DOWNSAMPLE; 
            cfg.detrend    = 'no';
            all_blocks     = ft_resampledata(cfg, all_blocks);
        end

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
        cfg1 = [];
        cfg1.length  = 4;
        cfg1.overlap = 0;
        all_blocks_segmented = ft_redefinetrial(cfg, all_blocks);
        %}

        cfg = [];
        cfg.output  = 'pow';
        cfg.channel = 'all';
        cfg.method  = 'mtmfft';
        cfg.taper   = 'boxcar';
        cfg.foi     = 0:0.01:30; % 1 / cfg1.length = 0.25 (the longer the segments, the more reso we can have here)
                                  % so for a reso of 0.01Hz, we need at least 1 segment with a length of 1 / 0.01 = 100 seconds
        freq         = ft_freqanalysis(cfg, all_blocks);

        freq.freq   = cfg.foi; % for some reason the "freq" field contains non-whole numbers, fix it manually

        
        % where to save the figures
        save_location = [SubjectFolder 'Figures\\' run_name '\\'];
        mkdir(save_location);
        
        % this fn takes care of all the plotting 
        % (power spectrum & topo for each freq band)
        plot_TFR(freq, lay, save_location, [1 30], true);
        
        
        % SAVE all relevant variables from the workspace
        save(S3_output_file, 'SubjectFolder', 'run_name', 'freq', '-v7.3');       
    else
        load(S3_output_file);
    end
end
