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
common_EM();

% location to save the results for all subjects
mkdir(ResultsFolder_thisrun);
mkdir(ResultsFolder_conn_thisrun);

% load layout & neighbours
load(LAYOUT_FILE); % use our custom-made layout & neighbours
%figure; ft_plot_layout(lay);
load(NEIGHBOURS_FILE);
load(ALL_LABELS_FILE); % list of real EEG channels (i.e. excluding EOG & ref chans): 61 chans for AntNeuro, 27 chans for NeuroPrax

% start EEGLAB if needed
if PLOT_CHANNEL_SPECTRA
    eeglab;
end


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
        - ICA artefact rejection (optional)
        - channel repair / interpolation
        - offline rereferencing
        % trigger-based trial definition (i.e. epoching)
        %}
        
        
        % get the eeg data file name 
        files = dir(fullfile(SubjectFolder, CONFILE_NAME));
        rawfile = fullfile(files(1).folder, files(1).name);

        
        % >>>
        % Step 1: filtering
        output_file = [output_path 'filtered.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2) 
            % first, read in the raw data
            % (keep this code here - rather than moving into filtering.m -
            % coz the cfg for ft_definetrial may need to be customised for each dataset)
            
            %hdr = ft_read_header(rawfile);%, 'dataformat','yokogawa_con'); % read header file         

            % ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
            cfg                      = [];
            cfg.trialfun             = 'ft_trialfun_general';
            cfg.datafile             = rawfile;
            %cfg.headerfile           = [rawfile(1:end-3) 'vhdr'];
            cfg.trialdef.triallength = Inf;
            cfg.trialdef.ntrials     = 1; % read in all data as a single segment, coz filtering should be done on continuous data
            cfg = ft_definetrial(cfg);

            % https://www.fieldtriptoolbox.org/workshop/madrid2019/tutorial_cleaning/
            cfg.demean     = 'yes';
            cfg.detrend    = 'no';
            cfg.continuous = 'yes';
            alldata = ft_preprocessing(cfg);                       
            
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
        
        %{
        % TROUBLESHOOTING of flat channels!
        % plot fft
        cfg = [];
        cfg.output  = 'pow';
        cfg.channel = 'all';
        cfg.method  = 'mtmfft';
        cfg.taper   = 'boxcar';
        cfg.foi     = 1:1:50; % 1 / cfg1.length = 0.25 (the longer the segments, the more reso we can have here)
                                  % so for a reso of 0.005Hz, we need at least 1 segment with a length of 1 / 0.005 = 200 seconds
        freq         = ft_freqanalysis(cfg, alldata);

        figure; hold on;
        channels = 1:27;
        %channels([6 7 14 22 26]) = [];
        %channels = [9 25];
        for chan = channels
            plot(freq.freq, freq.powspctrm(chan,:))
        end
        xlim([1 50]);
        xlabel('Frequency (Hz)');
        ylabel('Absolute power (uV^2)');
        hold off;
        %}
        
        
        % >>>
        % Step 2: do some additional processing (for migraine study only)        

        % For NeuroPrax EEG, remove the prefix "EEG " from channel labels
        % so that later on we can directly type in channel labels in the
        % ft_rejectvisual GUI (cannot enter directly if label contains space)
        if contains(alldata.label(1), 'EEG')
            temp = alldata.label(1:27);
            temp = cellfun(@(x) x(5:end), temp, 'un', 0); % remove first 4 chars in each cell
            alldata.label(1:27) = temp;
            alldata.hdr.label(1:27) = temp; % also update this (just in case)
        end
        
        % Extract the EC / EO section of the recording:
        % (only do this for old data - Bec now records EC & EO into 2
        % separate files, and each file is 15mins - we manually select the
        % best 5mins by marking the rest as artefact)
        if ~strcmp(SubjectID(1:9), 'Subject_2') % Subject_2xx would be Bec's data
            if strcmp(SUBJ_GROUP, 'patients') % for migraineurs, all subjects had the same timing
                if strcmp(run_name(1:3), '_EC')
                    latency = [0 290]; % EC: 0 - 290 sec
                elseif strcmp(run_name(1:3), '_EO')
                    latency = [300 10000]; % EO: from 300 sec onwards
                end
            elseif strcmp(SUBJ_GROUP, 'controls') % for controls, look up the timing we manually created and saved for each subject
                load([SubjectFolder 'EC_EO_timing.mat']);
                if strcmp(run_name(1:3), '_EC')
                    latency = EC_EO_timing.EC;
                elseif strcmp(run_name(1:3), '_EO')
                    latency = EC_EO_timing.EO;
                end
            end
            % now keep the relevant section of data & discard the rest
            cfg = [];
            cfg.latency = latency;
            alldata = ft_selectdata(cfg, alldata);
        end

        
        % >>>
        % Step 3: manually mark artefact sections
        output_file = [output_path 'arft.mat'];
        
        % if haven't already processed this before, do it now & save a copy
        if (exist(output_file, 'file') ~= 2)   
            % Print out SubjectID so we know which subject we are working on
            fprintf(['\nCURRENT SUBJECT: ' SubjectID '\n\n']); 
            
            % print out the auto-detected noisy channels (based on RMS) 
            % this can be used as ref during the visual inspection
            list_NoisyChans = '';
            counter_NoisyChans = 0;
            for chan = 1:length(all_labels) % only check the 27 real channels
                if rms(alldata.trial{1,1}(chan,:)) > 35 % threshold for EC: 35
                    list_NoisyChans = [list_NoisyChans ' ' int2str(chan)]; % add chan to list
                    counter_NoisyChans = counter_NoisyChans + 1; % increment counter
                end
            end
            fprintf('Suggest removing at least %d channels (rms > 35): %s\n', counter_NoisyChans, list_NoisyChans);

            % TROUBLESHOOTING - plot channel spectra (from raw data) using eeglab
            if PLOT_CHANNEL_SPECTRA
                EEG = pop_biosig(rawfile, 'blockrange',[0 290]); 
                EEG = eeg_checkset( EEG );
                % read chanlocs & remove non-EEG channels
                if length(EEG.chanlocs) == 33
                    EEG = pop_chanedit(EEG, 'load',{'Z:\\Analysis\\Judy\\EEG_analysis_scripts\\chanlocs_XYZ_33chan_forEEGLAB.txt' 'filetype' 'sfp'});
                elseif length(EEG.chanlocs) == 32
                    EEG = pop_chanedit(EEG, 'load',{'Z:\\Analysis\\Judy\\EEG_analysis_scripts\\chanlocs_XYZ_32chan_forEEGLAB.txt' 'filetype' 'sfp'});
                end
                EEG = pop_select( EEG,'nochannel',{'EOG1' 'EOG2' 'REF' 'ECG1' 'ECG2' 'TRIG'});
                EEG = eeg_checkset( EEG );
                % plot channel spectra
                figure('name',[SubjectID ' - spectopo()']); 
                pop_spectopo(EEG, 1, [0  289998], 'EEG' , 'freq', [10 30 50], 'freqrange',[2 60],'electrodes','off');
                EEG = eeg_checkset( EEG );
            end
            
            [arft] = mark_artefact(alldata, [-180 180]);
            
            if ~BROWSING_WITHOUT_SAVE
                save(output_file, 'arft', '-v7.3');
            end
        else
            load(output_file);
        end       
        
        % If running in batch, skip to next subject now
        if BROWSING_WITHOUT_SAVE
            continue;
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
        output_file = [output_path 'selChLabel' file_suffix '.mat'];
                
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
            
            if ~BROWSING_WITHOUT_SAVE
                save(output_file, 'selChLabel', '-v7.3');
            end
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
            output_file = [output_path 'afterICA' file_suffix '.mat'];
        
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
%
for i = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    
    output_path = [SubjectFolder output_name];
    S3_output_file = [ResultsFolder_thisrun SubjectID S3_output_filename];

    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S3_output_file, 'file') ~= 2)

        %load([output_path S2_output_filename]);
        load([output_path S1_output_filename]);
        
        % where to save the figures
        save_location = [SubjectFolder 'Figures' run_name '\\'];
        mkdir(save_location);
        % add prefix to filename if needed
        if ~isempty(file_suffix)
            save_location = [save_location file_suffix(2:end) '_'];
        end
        
        
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
        cfg.foi     = 1:1:30; % 1 / cfg1.length = 0.25 (the longer the segments, the more reso we can have here)
                                  % so for a reso of 0.005Hz, we need at least 1 segment with a length of 1 / 0.005 = 200 seconds
        freq        = ft_freqanalysis(cfg, all_blocks);
        
        freq.freq   = cfg.foi; % for some reason the "freq" field contains non-whole numbers, fix it manually

        
        % this fn takes care of all the plotting 
        % (power spectrum & topo for each freq band)
        plot_TFR(freq, lay, save_location, [1 30], FREQ_BANDS);
        
        
        % SAVE all relevant variables from the workspace
        save(S3_output_file, 'SubjectFolder', 'run_name', 'freq', '-v7.3');       
    else
        load(S3_output_file);
    end
end


%% Stage 4: connectivity analysis

for i = 1:length(SubjectIDs)
    
    SubjectID = cell2mat(SubjectIDs(i));
    SubjectFolder = [DataFolder SubjectID '\\'];
    
    output_path = [SubjectFolder output_name];
    S4_output_file = [ResultsFolder_conn_thisrun SubjectID S4_output_filename];

    % if haven't already processed this stage before, do it now & save a copy
    if (exist(S4_output_file, 'file') ~= 2)

        load([output_path S1_output_filename]);
 
        % where to save the figures
        save_location = [SubjectFolder 'Figures' run_name '\\connectivity\\'];
        mkdir(save_location);
        % add prefix to filename if needed
        if ~isempty(file_suffix)
            save_location = [save_location file_suffix(2:end) '_'];
        end

        
        % Apply surface Laplacian (SL) to deal with volumn conduction issue
        % https://www.youtube.com/watch?v=CodQ5-pmXdQ
        % https://www.fieldtriptoolbox.org/reference/ft_scalpcurrentdensity/
        if APPLY_SL
            cfg = [];
            cfg.method = 'spline'; %'finite'; % finite-difference method for the surface Laplacian on a triangulated sphere
            cfg.order = 4; % use less flexible splines (4 or 5) for low-density EEG (e.g. 32 chan). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4537804/
            %cfg.feedback = string, 'no', 'text', 'textbar', 'gui' (default = 'text')
            load(ELEC_FILE); % just load the version we have already made
            cfg.elec = elec;
            [all_blocks] = ft_scalpcurrentdensity(cfg, all_blocks);
        end
        
        
        % (1) Coherence
        % https://www.fieldtriptoolbox.org/tutorial/connectivity/#non-parametric-computation-of-the-cross-spectral-density-matrix
        cfg           = [];
        cfg.method    = 'mtmfft';
        cfg.taper     = 'dpss';
        cfg.output    = 'fourier';
        cfg.tapsmofrq = 2; % frequency smoothing of 2Hz
        cfg.foi       = 1:1:30; 
        %cfg.pad       = 'nextpow2'; % for more efficient FFT computation (but sometimes worse)
        freq          = ft_freqanalysis(cfg, all_blocks);
        
        cfg           = [];
        cfg.method    = 'coh';
        coh           = ft_connectivityanalysis(cfg, freq);
        %cohm          = ft_connectivityanalysis(cfg, mfreq);
        
        cfg           = [];
        cfg.parameter = 'cohspctrm';
        cfg.xlim      = [2 30]; % we are interested in 2-30Hz (everything else was filtered out)
        cfg.zlim      = [0 1];
        ft_connectivityplot(cfg, coh); %, cohm);

        set(gcf, 'Position', get(0, 'Screensize')); % make the figure full-screen
        if APPLY_SL
            export_fig(gcf, [save_location 'coherence_afterSplineSL.png']); % use this tool to save the figure exactly as shown on screen
        else
            export_fig(gcf, [save_location 'coherence.png']); % use this tool to save the figure exactly as shown on screen
        end
        
        
        % (2) Granger causality (directional connectivity)
        % https://www.fieldtriptoolbox.org/tutorial/connectivity/#computation-of-the-multivariate-autoregressive-model
        
        % Only works if data are continuous (or epochs of same length?),
        % so any subjects with manually marked bad sections need to be concatenated first
        % Tried running on subject 661 - see figure in its connectivity folder
        %{
        % compute the multivariate autoregressive model
        cfg         = [];
        cfg.order   = 2; %5;
        cfg.toolbox = 'biosig'; % 'bsmart' toolbox didn't work
        mdata       = ft_mvaranalysis(cfg, all_blocks);
        
        % compute the spectral transfer function
        cfg        = [];
        cfg.method = 'mvar';
        mfreq      = ft_freqanalysis(cfg, mdata);
        
        % compute & plot the granger causality spectrum
        cfg           = [];
        cfg.method    = 'granger';
        granger       = ft_connectivityanalysis(cfg, mfreq);

        cfg           = [];
        cfg.parameter = 'grangerspctrm';
        cfg.xlim      = [2 30]; % we are interested in 2-30Hz (everything else was filtered out)
        %cfg.zlim      = [0 1];
        ft_connectivityplot(cfg, granger);

        set(gcf, 'Position', get(0, 'Screensize')); % make the figure full-screen
        export_fig(gcf, [save_location 'granger.png']); % use this tool to save the figure exactly as shown on screen
        %}
        
        
        % SAVE all relevant variables from the workspace
        save(S4_output_file, 'SubjectFolder', 'run_name', 'freq', 'coh', '-v7.3');       
    else
        load(S4_output_file);
    end
end
