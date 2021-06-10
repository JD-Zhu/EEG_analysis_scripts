function source_v1
    % NOTE: always run from the "analysis_scripts" folder, to ensure all relative
    % paths work as expected
    cd('E:\Judy\Exp2\7_MEG-analysis\scripts');

    addpath(pwd); % allow access to scripts in this folder even after we cd into each SubjectFolder
    %addpath([pwd '\\coreg-master\\']); % allow access to coreg scripts

    % run the #define section
    global DataFolder; global ResultsFolder; global ResultsFolder_ROI; global ResultsFolder_Source;
    global filename_suffix; 
    global eventnames_real;
    common();

    % specify all paths as absolute paths here, because below we 'cd' into each
    % subject folder & the relative path will no longer be correct
    templates_dir = [pwd '\\..\\FT_templates\\']; % Location of the required templates (headmodel, grid & mri).
                                           % These usually come with the FT toolbox, but for some reason 
                                           % I don't have them (later found out: I had the "lite" version of FT).
                                           % For ease of access (consistent path across computers), I've stored a copy here.                                  

    % location of your MRI database (consistent relative path across computers)
    MRI_path = [pwd '\\..\\..\\..\\MRI_databases\\'];
    
    
    %% = Settings =
    % Please adjust as required:

    % run the source localisation (of sensor effects)?
    RUN_SOURCE_LOCALISATION = false;
    
    % run the ROI analysis?
    RUN_ROI_ANALYSIS = true;
    
    % use fixed or free dipole orientation in the beamformer?
    FIXED_ORI = 'no'; % 'yes' == fixed, 'no' = free
    
    % use which method ('svd' or 'centroid') for collapsing 
    % all vertices in an ROI into a single virtual sensor?
    % Note: if you set FIXED_ORI = 'no' above, a diff VE methods has to be
    % used (i.e. this var will not be looked at)
    VE_METHOD = 'svd';
    
    % which version of MEMES to use?
    % (MEMES1 - for comparability with earlier results; MEMES3 - for Chinese MRI database)
    MEMES_VERSION = 'MEMES3';

    
    % which preprocessing/ERF results to use?
    run_name = 'TSPCA10000_3';
    if strcmp(FIXED_ORI, 'yes')
        run_suffix = '';
    else % free dipole orientation
        run_suffix = '_freeori';    
    end
    ResultsFolder_thisrun = [ResultsFolder run_name '\\']; % ERF results for all subjects
    ResultsFolder_ROI_thisrun = [ResultsFolder_ROI run_name run_suffix '\\']; % ERF results for all subjects
    ResultsFolder_Source_thisrun = [ResultsFolder_Source run_name run_suffix '\\']; % ERF results for all subjects

    
    % = File names for saving =
    
    % filename for saving the beamformer output (to avoid running the whole thing every time)
    if strcmp(FIXED_ORI, 'yes') % fixed dipole orientation
        %Beamformer_output_filename = 'beamformer_Chinese.mat'; % This is an old version, where cfg.lcmv.lamda = '5%'; 
        Beamformer_output_filename = 'beamformer_Chinese_lambda=1.mat'; 
    else % free dipole orientation
        Beamformer_output_filename = 'beamformer_freeori.mat'; 
    end
    % Note: To save a different version of beamformer results (e.g. when using 
    % a new set of ERF outputs), simply change this filename.
    % Similarly, to load a previous version of beamformer results,
    % change this filename accordingly.
    
    % filename for saving the ROI activity (stored in ResultsFolder_ROI for all subjects)
    ROI_output_filename = '_ROI.mat';
    
    
    %% check input & settings are valid
    if ~(strcmp(FIXED_ORI, 'yes') || strcmp(FIXED_ORI, 'no'))
        fprintf('Error in source_v1: Invalid setting for FIXED_ORI\n');
    end
    

    %% start
    
    %SubjectFolders = listFolders(DataFolder);
    SubjectFolders = [dir([DataFolder 'A*']); dir([DataFolder 'B*'])];
    SubjectFolders = {SubjectFolders.name}; % extract the names into a cell array
    
    % each cycle processes one subject
    for h = 1:length(SubjectFolders)

        SubjectID = SubjectFolders{h};
        SubjectFolder = [DataFolder, '\\', SubjectID];
        cd(SubjectFolder); % change directory

        %% Step 1: coreg (to obtain individualised headmodel & sourcemodel)
        
        % check whether this is a subject folder
        elp_file  = dir('*.elp'); % find the .elp file
        if isempty(elp_file) % not a subject folder, skip
            fprintf('\n"%s": Not a subject folder - skip!\n', SubjectID);
            continue;
        end

        % prepare all the inputs to MEMES
        filename_base = elp_file.name(1:strfind(elp_file.name,'.')-1); % get the base filename (ie. remove suffix)
        elpfile = [filename_base, '.elp'];
        hspfile = [filename_base, '.hsp'];
        confile = [filename_base, '_B1.con'];
        mrkfile = [filename_base, '_ini.mrk']; % choose which marker file to use
        
        % specify the bad marker coils (max 2) for each subject
        % Enter as: {'LPAred','RPAyel','PFblue','LPFwh','RPFblack'}
        bad_coil = ''; 
        if strcmp(SubjectID, 'A02-EL-3604') || strcmp(SubjectID, 'A05-RW-3584') ...
          || strcmp(SubjectID, 'A11-DZ-3541') || strcmp(SubjectID, 'B07-OS-3547') ...
          || strcmp(SubjectID, 'B02-YW-3523') % B02 & B07 had 3 markers that moved >5mm, we can only exclude 2 here (however, the movement was only btwn preB2 & post, so it could have happened at the very end after task finished)
            bad_coil = {'LPAred', 'PFblue'};
        elseif strcmp(SubjectID, 'A03-ZZ-3555') || strcmp(SubjectID, 'B03-BQ-3554')
            bad_coil = {'LPAred'};
        end

                
        % check which version of MEMES to use
        if strcmp(MEMES_VERSION, 'MEMES1')
            MRI_folder = [MRI_path 'HCP\\']; % HCP database
            coreg_output = [pwd '\\MEMES\\']; % where to store the output from MEMES
            
            % if headmodel etc haven't been generated, do this now
            if ~exist([coreg_output 'headmodel_singleshell.mat'], 'file')
                % prepare extra inputs to MEMES
                temp = load([MRI_folder 'mesh_library.mat']);
                mesh_library = temp.mesh_library;
                initial_mri_realign = temp.initial_mri_realign;
                path_to_MRI_library = MRI_folder;

                MEMES(pwd,coreg_output,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,...
                    mesh_library,initial_mri_realign,bad_coil);
                %MEMES2(pwd, elpfile, hspfile, confile, mrkfile,path_to_MRI_library,...
                %    mesh_library,initial_mri_realign, bad_coil, 'best', 1)

                %mrifile = 'single_subj_T1.nii'; % use the template that comes with FT
                % don't use the dummy .mri file created by ME160 - it doesn't contain anything
                %coreg_yokogawa_icp(pwd, confile, mrkfile, mrifile, ...
                %    hspfile, elpfile, 100, 0.1); %[pwd,'\\',filename_base,'_MRI\\',filename_base,'.nii'], ...
            end
            
            % now load the headmodel & grads & mri generated by MEMES
            temp = load ([coreg_output 'headmodel_singleshell.mat']);
            headmodel = temp.headmodel_singleshell;
            temp = load ([coreg_output 'grad_trans.mat']);
            grads = temp.grad_trans;
            
            if ~exist([coreg_output 'mri_realigned_transformed.mat'], 'file')
                temp = load ([coreg_output 'mri_realigned.mat']);
                temp1 = load ([coreg_output 'trans_matrix.mat']);
                mri_realigned = ft_transform_geometry(temp1.trans_matrix, temp.mri_realigned);
                save ([coreg_output 'mri_realigned_transformed.mat'], 'mri_realigned');
                mri = mri_realigned;
            else
                temp = load ([coreg_output 'mri_realigned_transformed.mat']);
                mri = temp.mri_realigned;
            end
            
        elseif strcmp(MEMES_VERSION, 'MEMES3')
            MRI_folder = [MRI_path 'SLIM_completed\\']; % SLIM Chinese database
            %MRI_folder = [MRI_path 'HCP_for_MEMES3\\']; % new HCP database (works with MEMES3)
            coreg_output = [pwd '\\MEMES_best\\']; % where to store the output from MEMES
            %coreg_output = [pwd '\\MEMES3_Chinese\\'];
            %coreg_output = [pwd '\\MEMES3_HCP\\'];        
            
            % if headmodel etc haven't been generated, do this now
            if ~exist([coreg_output 'headmodel.mat'], 'file')
                
                % call MEMES3 (new version Oct 2019)
                if length(bad_coil) >= 2     % we use 'rot3dfit' by default (faster), however this sometimes has issues (creates upside-down coreg)
                    realign_method = 'icp';  % when there are 2 or more bad coils, so in that case we use 'icp' instead
                else
                    realign_method = 'rot3dfit';
                end
                [headshape_downsampled] = downsample_headshape(hspfile, 'no', 0); % Robert recommends this setting for my defaced MRI database
                [grad_trans] = mq_realign_sens(pwd,elpfile,hspfile,confile,mrkfile,bad_coil,realign_method);
                MEMES3(pwd, grad_trans, headshape_downsampled, MRI_folder, 'best', [0.98:0.01:1.12], 5, []); % do not set a weighting for facial information
                
                % call MEMES3 (old version 2018)
                %MEMES3_old_2018(pwd, elpfile, hspfile, confile, mrkfile, MRI_folder, bad_coil, 'best', [0.99:0.01:1.01], 5, 'yes')
                
                
                % close the figures MEMES created (each subject creates 5
                % figures - becomes too many when running in batch)
                close all;
                
                % move the MEMES output into the coreg_output folder.
                if ~exist(coreg_output, 'dir')
                    mkdir(coreg_output);
                end                
                % the files to move:
                % grad_trans, headshape, headmodel, trans_matrix, sourcemodel3d, shape
                movefile('*trans*', coreg_output);
                movefile('*shape*', coreg_output);
                movefile('*model*', coreg_output);
                movefile('*quality*', coreg_output);
                movefile('*example*', coreg_output);
                movefile('*scaling*', coreg_output);
                movefile('*realigned*', coreg_output);
                movefile('*winner*', coreg_output);
            end
 

            % now load the headmodel & grads & mri generated by MEMES
            temp = load ([coreg_output 'headmodel.mat']);
            headmodel = temp.headmodel;
            temp = load ([coreg_output 'grad_trans.mat']);
            grads = temp.grad_trans;
            temp = load ([coreg_output 'sourcemodel3d.mat']);
            sourcemodel = temp.sourcemodel3d;

        else
            fprintf('\nError in source_v1.m: Incorrect selection of MEMES version. Please specify at the top of the script.\n');
        end   


        %% Step 2: load this subject's ERF results
        %
        subject_data = load([ResultsFolder_thisrun SubjectID '_erf' filename_suffix '.mat']);
        erf = subject_data.erf_clean;
        erf_cue_combined = subject_data.erf_allconds;
        %erf_target_combined = subject_data.erf_target_combined;
        clear subject_data;
        %
        
        %% only need this section if restricting to a certain freq band
        % (remember: erf is averaging then calc, freq is calc'ing on indi trials then average)
        % NOTE: you can do this elsewhere & simply load the results (just like you did above with the ERF results) 
    %{

        load trials_clean   % trials_clean is after visual rejection & just before ft_timelockanalysis in Yokogawa_FT_v5
        data = trials_clean.cue;

        % we saved the bad trials as NaN, so need to exclude these
        good_trials_idx = find(~isnan(cell2mat(cellfun(@(isgood)isgood(1),data.trial,'uni',0)))); %just need to evaluate the first element as all samples in bad trial are NaN

        cfg        = [];
        cfg.trials = good_trials_idx;
        data       = ft_redefinetrial(cfg, data);


        %% Bandpass Filter (to get the required freq band)
        cfg           = [];
        cfg.channel   = 'all';
        cfg.bpfilter  = 'yes';
        cfg.bpfreq    = [7 12];    % keep alpha only
        data_filtered = ft_preprocessing(cfg,data);


        %% Here we redefine trials based on the time-points of interest.
        % Make sure the timepoints are of equivalent length - no need to do this
        % but I'm leaving it here anyway
        cfg        = [];
        cfg.toilim = [-0.4 0]; % baseline window
        datapre    = ft_redefinetrial(cfg, data_filtered);

        cfg.toilim = [0.3 0.7]; % window of interest
        datapost   = ft_redefinetrial(cfg, data_filtered);

        % Here we are keeping all parts of the trial for your covariance matrix
        % **SHOULD THIS GO AFTER AVERAGE?***
        cfg                  = [];
        cfg.covariance       = 'yes';
        cfg.covariancewindow = [-1.0 1.0];
        avg                  = ft_timelockanalysis(cfg,data_filtered);

        % Time lock analysis for datapre and datapost period
        cfg                  = [];
        cfg.covariance       = 'yes';
        %cfg.covariancewindow = [-1.0 1.0];
        avgpre               = ft_timelockanalysis(cfg,datapre);
        avgpst               = ft_timelockanalysis(cfg,datapost);
    %}    


        %% if haven't already done the beamforming before, do it now & save a copy
        Beamformer_output_file = [pwd '\\' Beamformer_output_filename];
        if (exist(Beamformer_output_file, 'file') ~= 2)    

            %% Step 3: prepare sourcemodel & leadfield

            % check which version of MEMES we used
            % if MEMES1, then we need to create the individual sourcemodel here
            if strcmp(MEMES_VERSION, 'MEMES1')

                % Create template grid (aka. template sourcemodel)
                % Method 1:
                %{
                % Load template headmodel
                load([templates_dir, 'standard_singleshell']);
                template_headmodel = vol; % rename the loaded variable                                     
                template_headmodel = ft_convert_units(template_headmodel, 'cm'); % convert headmodel into standard units (CTF convention is to express everything in cm)

                % create template sourcemodel using template headmodel
                cfg             = [];
                cfg.grid.xgrid  = -20:0.5:20; % x-axis range: -20 ~ 20cm, step 0.5cm (change this value to adjust reso)
                cfg.grid.ygrid  = -20:0.5:20;
                cfg.grid.zgrid  = -20:0.5:20;
                cfg.grid.unit   = 'cm';
                cfg.grid.tight  = 'yes';
                cfg.inwardshift = -0.1;
                cfg.headmodel   = template_headmodel;
                template_sourcemodel = ft_prepare_sourcemodel(cfg);
                %}
                % Method 2: load the template sourcemodel provided by FT
                temp = load([templates_dir, 'standard_sourcemodel3d5mm']); % usually 10mm grid is fine.
                                                                    % here we use a higher resolution grid (5mm),
                                                                    % to increase the number of vertices assinged to
                                                                    % the parcel 'Frontal_Med_Orb_L' (in AAL atlas)
                                                                    % Otherwise it only contains 3 vertices - not enough to compute centroid!
                template_sourcemodel = temp.sourcemodel; % rename the loaded variable

                % Warp template grid into subject space
                cfg                = [];
                cfg.grid.warpmni   = 'yes';
                cfg.grid.template  = template_sourcemodel; % standard sourcemodel
                cfg.grid.nonlinear = 'yes';
                cfg.mri            = mri; % individual mri
                sourcemodel = ft_prepare_sourcemodel(cfg); % creates individual sourcemodel
                                                           % (the grid points map 1-to-1 onto the template grid points, with the .pos field 
                                                           % specifying the actual coordinates of these grid points in subject space)
            end                                 
            
            % if we used "MEMES2" or "MEMES3", then the individual grid
            % (incl. warping) has already been done for you. 
            % We loaded the desired grid (e.g. "sourcemodel3d_5mm") 
            % from the MEMES output above, so just plug it directly into ft_prepare_leadfield

            
            % Plots for sanity checks
            % plot headmodel, grid, and mri
            %{       
            ft_determine_coordsys(mri, 'interactive','no'); hold on
            ft_plot_vol(headmodel);
            ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
            %}           
            % plot headmodel, grid, and sensor locations
            %
            figure;
            ft_plot_sens(grads, 'style', '*b'); % plot the MEG sensor locations
            ft_plot_vol(headmodel, 'edgecolor', 'cortex'); alpha 0.4; % plot the single shell (i.e. brain shape)
            ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); % plot all vertices (ie. grid points) that are inside the brain
            %
            
            
            % Create the leadfield
            cfg            = [];
            cfg.grad       = grads;
            cfg.headmodel  = headmodel; % individual headmodel (from coreg)
            cfg.reducerank = 2; % Should check this is appropriate - also check the rank of the data as we project out mouth artifacts earlier
            cfg.channel    = erf_cue_combined.label; % use the actual channels present in our data (i.e. ensure that rejected sensors are also removed here)
            cfg.grid       = sourcemodel; % individual sourcemodel (warped from template grid)
            [grid]    = ft_prepare_leadfield(cfg); % sourcemodel + leadfield


            %% Step 4: LCMV beamformer
            % http://www.fieldtriptoolbox.org/tutorial/salzburg
            % http://www.fieldtriptoolbox.org/example/common_filters_in_beamforming

            % Create common spatial filter (based on data across all conds),
            % otherwise it would be circular (diff filters for diff conds, created based on each cond)
        %{    
            % first, compute the average erf (across 4 conds) for cue window & target window
            erf_cue_combined = erf.cuechstay; % average erf for cue window
            erf_cue_combined.avg = (erf.cuechstay.avg + erf.cuechswitch.avg + erf.cueenstay.avg + erf.cueenswitch.avg) / 4;
            erf_target_combined = erf.targetchstay; % average erf for target window
            erf_target_combined.avg = (erf.targetchstay.avg + erf.targetchswitch.avg + erf.targetenstay.avg + erf.targetenswitch.avg) / 4;
        %}    
            % run ft_sourceanalysis on the cue_combined erf & target_combined erf
            % to create common spatial filters (1 filter for cue, 1 filter for target)
            cfg                 = [];
            cfg.keeptrials      = 'no';
            cfg.channel         = 'MEG';
            cfg.grad            = grads;
            cfg.senstype        = 'MEG';
            cfg.method          = 'lcmv';
            cfg.grid            = grid; % individual sourcemodel + leadfield (warped from template grid)
            cfg.grid.unit       = 'cm';
            cfg.headmodel       = headmodel; % individual headmodel (from coreg)
            %cfg.lcmv.lamda      = '5%';
            cfg.lcmv.lamda      = '100%'; % regularisation parameter - set to 1 here coz our data is rank-reduced (due to rejecting ICA comps)
            cfg.lcmv.fixedori   = FIXED_ORI; % use the setting specified at the top
            cfg.lcmv.keepfilter = 'yes';
            cfg.lcmv.projectmom = 'no';
            cfg.lcmv.normalize  = 'yes'; %corrects for depth bias?
            source_cue_combined = ft_sourceanalysis(cfg, erf_cue_combined); % create spatial filter for cue window
            %source_target_combined = ft_sourceanalysis(cfg, erf_target_combined); % create spatial filter for target window

            % save
            save(Beamformer_output_file, 'sourcemodel', 'source_cue_combined'); % 'source_target_combined');
        end
        
                    
        %% Step 5: Source localisation (of sensor-space effects)
        % Note: ft_sourceanalysis is only for localisation (no time dimension, only creates a source image showing where in the 3d space is one cond different from another cond)
        
        if RUN_SOURCE_LOCALISATION
            % load required files
            temp = load([templates_dir, 'standard_sourcemodel3d5mm']);
            template_sourcemodel = temp.sourcemodel;

            temp = load(Beamformer_output_file);
            source_cue_combined = temp.source_cue_combined;
            source_target_combined = temp.source_target_combined;
            
            
            % time window to avg over (determined from sensor-space analysis)
            cue_ttype_window = [0.435 0.535];
            target_lang_window = [0.255 0.305];


            % Main effect of switch in cue window
            cue_stay = erf.cuechstay;
            cue_switch = erf.cuechswitch;
            cue_stay.avg = (erf.cuechstay.avg + erf.cueenstay.avg) / 2;
            cue_switch.avg = (erf.cuechswitch.avg + erf.cueenswitch.avg) / 2;

            cfg_select = []; % avgovertime in the time window of sensor effect
            cfg_select.latency     = cue_ttype_window;
            cfg_select.avgovertime = 'yes';
            cfg_select.nanmean     = 'yes';
            cue_stay = ft_selectdata(cfg_select, cue_stay);
            cue_switch = ft_selectdata(cfg_select, cue_switch);

            cfg.grid.filter = source_cue_combined.avg.filter; % use the common spatial filter created above
            source_cue_stay = ft_sourceanalysis(cfg, cue_stay);
            source_cue_switch = ft_sourceanalysis(cfg, cue_switch);

            % Main effect of lang in target window
            target_ch = erf.targetchstay;
            target_en = erf.targetenstay;
            target_ch.avg = (erf.targetchstay.avg + erf.targetchswitch.avg) / 2;
            target_en.avg = (erf.targetenstay.avg + erf.targetenswitch.avg) / 2;

            cfg_select = []; % avgovertime in the time window of sensor effect
            cfg_select.latency     = target_lang_window;
            cfg_select.avgovertime = 'yes';
            cfg_select.nanmean     = 'yes';
            target_ch = ft_selectdata(cfg_select, target_ch);
            target_en = ft_selectdata(cfg_select, target_en);

            cfg.grid.filter = source_target_combined.avg.filter; % use the common spatial filter created above
            source_target_ch = ft_sourceanalysis(cfg, target_ch);
            source_target_en = ft_sourceanalysis(cfg, target_en);


            % load template MRI for warping into common (MNI) space
            template_mri          = ft_read_mri(fullfile(templates_dir, 'single_subj_T1_1mm.nii')); %This is the standard which matches AAL space
            template_mri.coordsys = 'nifti_spm'; % so that FieldTrip knows how to interpret the coordinate system

            % localise the source of each effect (ie. location in the brain that shows strongest diff btwn the 2 conds)
            save_filename = [ResultsFolder_Source_thisrun 'cue_ttype\\' SubjectID];
            localise_effect_source(source_cue_stay, source_cue_switch, template_sourcemodel, template_mri, save_filename);

            save_filename = [ResultsFolder_Source_thisrun 'target_lang\\' SubjectID];
            localise_effect_source(source_target_ch, source_target_en, template_sourcemodel, template_mri, save_filename);
        end


        %% Step 6: ROI analysis (source reconstruction) (independent from sensor-level results)
        % Here we use the atlas to create VEs (virtual sensors) - 1 VE represents 1 ROI

        if RUN_ROI_ANALYSIS
            
            % if haven't done ROI activity reconstruction before, do it now & save a copy
            ROI_output_file = [ResultsFolder_ROI_thisrun SubjectID ROI_output_filename];
            if exist(ROI_output_file, 'file') ~= 2
                % load required files
                temp = load([templates_dir, 'standard_sourcemodel3d5mm']);
                template_sourcemodel = temp.sourcemodel;
                template_sourcemodel = ft_convert_units(template_sourcemodel, 'mm');

                temp = load(Beamformer_output_file);
                sourcemodel = temp.sourcemodel;
                source_cue_combined = temp.source_cue_combined;
                %source_target_combined = temp.source_target_combined;


                % Load Atlas (contains parcellation of brain into regions/tissues/parcels)
                atlas = ft_read_atlas(fullfile(templates_dir, 'ROI_MNI_V4.nii'));
                atlas = ft_convert_units(atlas, 'mm');% ensure that atlas and template_sourcemodel are expressed in the same units

                % Interpolate the atlas onto template sourcemodel (10mm grid),
                % because the atlas may not be at the same resolution as your grid
                % (e.g. you created a grid with 6000 vertices, but atlas may only have 2000 vertices)
                cfg                  = [];
                cfg.interpmethod     = 'nearest';
                cfg.parameter        = 'tissue';
                atlas_interpo = ft_sourceinterpolate(cfg, atlas, template_sourcemodel);

                % Add the tissue labels from atlas to our sourcemodel
                % (each "tissue label" defines one "parcel")
                %sourcemodel.tissue      = atlas_interpo.tissue; 
                %sourcemodel.tissuelabel = atlas_interpo.tissuelabel;


                % Define our ROIs (can combine multiple parcels together to form one ROI)
                ROIs = {{'Frontal_Inf_Oper_L';'Frontal_Inf_Tri_L'},{'Frontal_Inf_Oper_R';'Frontal_Inf_Tri_R'},...
                    {'Temporal_Sup_L'},{'Temporal_Sup_R'},{'Supp_Motor_Area_L'},{'Supp_Motor_Area_R'},...
                    {'Cingulum_Ant_L'},{'Cingulum_Ant_R'},{'Frontal_Mid_L'},{'Frontal_Mid_R'},...
                    {'SupraMarginal_L'},{'SupraMarginal_R'},{'Calcarine_L';'Calcarine_R'}};
               
                ROIs_label = {'LIFG','RIFG','LSTG','RSTG','LSMA','RSMA','LACC','RACC','LdlPFC','RdlPFC','LSMG','RSMG','V1'}; %Labels for the groupings
                % For dlPFC (BA9, 10, 46) <- I chose 'Frontal_Mid_L' (middle frontal gyrus)
                % According to Wikipedia, MFG == BA9,10,46 
                % (https://en.wikipedia.org/wiki/Middle_frontal_gyrus)
                % Also, "The DLPFC is not an anatomical structure, but rather a functional one. 
                % It lies in the middle frontal gyrus of humans." 
                % (https://en.wikipedia.org/wiki/Dorsolateral_prefrontal_cortex)
                %
                % Have tried 'Frontal_Med_Orb_L' (medial orbitofrontal cortex) before, 
                % the location is the medial-orbital side of the frontal tip of IFG.
                % 
                % Some advice re: ROI definitions
                % (https://sourceforge.net/p/marsbar/mailman/message/34517456/)
                % "I have seen many students trying to replicate prior work on “DLPFC” and other nonsense terms, and usually fail, because there usually aren't any protocols - or consensus.
                % I suggest you do stick to anatomical terms. Another possibility may be to define your ROI functionally (which you can do easily in MarsBaR, starting from SPM).
                % Incidentally, AAL is a single-subject atlas. Any single subject atlas is on average about five points less accurate compared with the manual gold standard on standard overlap indices (e.g. Dice) than a maximum probability atlas based on multiple subjects, simply because any individual is not representative... Maximum probability atlases are single files, hence exactly as easy to handle as AAL, and readily available, e.g. from LONI, or from us. Available at http://www.brain-development.org "
                % 
                %
                % For more precise definition of preSMA, can try:
                % (1) Stanford FIND parcellation (2) HMAT
                %

        
                % Make a plot showing the vertices in the parcels on the source model - a good sanity check 
                %{
                % to plot selected ROIs only:
                %ROIs = {{'Frontal_Med_Orb_L'},{'Frontal_Med_Orb_R'},{'Frontal_Mid_L'},{'Frontal_Mid_R'}};       

                % draw the ROIs
                parcel_vertices = [];
                for e = 1:length(ROIs)
                    % make a separate figure for each ROI
                    figure('Name','Position of Points','NumberTitle','off'); hold on;
                    ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
                    hold on;

                    indx_pos = [];
                    for region = find(ismember(atlas_interpo.tissuelabel, char(ROIs{e})))
                        % Get atlas points
                        indx = (find(atlas_interpo.tissue == region));
                        indx_pos_temp = sourcemodel.pos(indx,:);
                        %for vol = 1:length(indx)     % this is incorrect
                        %    indx_pos_temp(vol,:) = sourcemodel.pos(indx(vol),:); % take the warped coords
                        %end
                        indx_pos = vertcat(indx_pos, indx_pos_temp);
                    end
                    ft_plot_mesh(indx_pos, 'vertexcolor','k', 'vertexsize',10); hold on;
                    title(atlas_interpo.tissuelabel(region), 'Interpreter', 'none'); 
                    
                    % create a struct with the vertices for each ROI being used. Might be useful at some point
                    parcel_vertices.(ROIs_label{e}) = indx_pos; 
                end
                %}


                %% Create VE for each ROI

                % Each cycle deals with one ROI
                for k = 1:length(ROIs)
                    ROI_name = ROIs_label{k};

                    % for this ROI, find a list of vertices that belong to it, and
                    % extract the spatial filter for each vertex in cue window & target window
                    vertices_all = []; % will hold a single list of all vertices (from all parcels belonging to this ROI)
                    for j = 1:length(ROIs{k})
                        indx        = find(ismember(atlas_interpo.tissuelabel, ROIs{k}{j})); % find index of the required tissue label
                        vertices    = find(atlas_interpo.tissue == indx); % find vertices that belong to this tissue label
                        % add vertices from the current parcel to the overall list
                        vertices_all = [vertices_all; vertices];
                    end
                    % get the spatial filter (i.e. set of weights) for each vertex
                    if strcmp(FIXED_ORI, 'yes') % for fixed orientation, there is 1 set of weights for each vertex
                        vertices_filters_cue = cat(1, source_cue_combined.avg.filter{vertices_all}); 
                    else % for free orientation, there are 3 sets of weights (1 on each axis) for each vertex
                        vertices_filters_cue = cat(3, source_cue_combined.avg.filter{vertices_all}); 
                    end

                    % create virtual sensor for this ROI, using the appropriate fn 
                    % according to what settings we selected at the top
                    if strcmp(FIXED_ORI, 'no') % free dipole orientation
                        VE = create_virtual_sensor_freeori(ROI_name, vertices_filters_cue, erf_cue_combined, erf, 1:length(eventnames_real)); 
                    elseif strcmp(VE_METHOD, 'centroid')
                        VE = create_virtual_sensor_Centroid(ROI_name, vertices_all, vertices_filters_cue, erf_cue_combined, erf, 1:length(eventnames_real), headmodel, sourcemodel);
                    else
                        VE = create_virtual_sensor_SVD(ROI_name, vertices_filters_cue, erf_cue_combined, erf, 1:length(eventnames_real)); 
                    end

                    if ~isempty(VE) % successful
                        ROI_activity.(ROI_name) = VE;
                    else
                        fprintf(['No solution for ', ROI_name, ' in cue window.']);
                    end

                    %{
                    % create virtual sensor for this ROI in target window
                    if (strcmp(VE_METHOD, 'centroid'))
                        VE = create_virtual_sensor_Centroid(ROI_name, vertices_all, vertices_filters_target, erf_target_combined, erf, conds_target, headmodel, sourcemodel);
                    else
                        VE = create_virtual_sensor_SVD(ROI_name, vertices_filters_target, erf_target_combined, erf, conds_target);
                    end

                    if ~isempty(VE) % successful
                        for j = conds_target  % append to existing cue-window results
                            ROI_activity.(ROI_name).(eventnames_real{j}) = VE.(eventnames_real{j});
                        end
                    else
                        fprintf(['No solution for ', ROI_name, ' in target window.']);
                    end
                    %}
                end

                save(ROI_output_file, 'ROI_activity');

                % Plot the source activity at each ROI
                %{
                for k = 1:length(ROIs)
                    ROI_name = (ROIs_label{k});

                    % plot for cue window
                    figure; hold on; 
                    title(['Cue window: ' ROI_name]);
                    for j = conds_cue
                       plot(ROI_activity.(ROI_name).(eventnames_real{j}).time, ROI_activity.(ROI_name).(eventnames_real{j}).avg);
                       xlim([-0.2 0.75]); % epoch was [-1 1], we only want to plot [-0.2 0.75]
                    end
                    legend(eventnames(conds_cue));

                    % plot for target window
                    figure; hold on; 
                    title(['Target window: ' ROI_name]);
                    for j = conds_target
                       plot(ROI_activity.(ROI_name).(eventnames_real{j}).time, ROI_activity.(ROI_name).(eventnames_real{j}).avg);
                       xlim([-0.2 0.75]); % epoch was [-1 1], we only want to plot [-0.2 0.75]
                    end
                    legend(eventnames(conds_target));
                end
                %}
            end

            % Statistical analysis on the ROI activities will be carried out in stats_ROI.m
        end
    end

    
    %% Grand average for Step 5 (Source localisation)
    % Q: do we actually do any stats here?
    % A: no, because we have already done the stats at sensor level, now we
    % just want to know where (in source space) the already-discovered effect occurs
    
    if RUN_SOURCE_LOCALISATION
    
        % get the list of contrasts we are examining (based on saved folder structure)
        contrasts = listFolders(ResultsFolder_Source);

        % each cycle processes one contrast (e.g. cue_ttype, target_lang)
        for index = 1:length(contrasts)
            % read in the saved blob for each subject (all in common space)
            blobs_folder = [ResultsFolder_Source cell2mat(contrasts(index)) '\\'];
            blobs_files = dir([blobs_folder 'M*.mat']);        
            for subject = 1:length(blobs_files)
                temp = load([blobs_folder blobs_files(subject).name]);
                blobs(subject) = temp.source_int;
                %blobs{subject} = temp.sourceDiff; % use this option if using z-scores
            end

            % average the blob across all subjects
            grandave = blobs(1);
            %grandave.pow = median([blobs.pow], 2); % use median   
            %save_filename = [blobs_folder 'average_blob--median'];
            grandave.pow = trimmean([blobs.pow], 12.5, 2); % use robust mean (removing top 1 & bottom 1 subject)      
            save_filename = [blobs_folder 'average_blob--trimmean'];

            % only retain vertices in the top 25%, set the rest to 0 (only need this step if the intensity threshold in xjview GUI doesn't seem to correspond with the actual intensity scale)
            threshold = quantile(grandave.pow, 0.75);
            grandave.pow(grandave.pow < threshold) = 0;

            % average the blob (z-scores) across all subjects
            %{
            cfg = [];
            cfg.parameter = 'zscores';
            [grandave] = ft_sourcegrandaverage(cfg, blobs{:});
            %}

            grandave = rmfield(grandave, 'cfg'); % useless field taking >2Gb space
            save(save_filename, 'grandave', '-v7.3');

            % smear onto mri overlay (can be diff resolution)
            %{
            cfg_source              = [];
            cfg_source.voxelcoord   = 'no';
            cfg_source.parameter    = 'zscores';
            cfg_source.interpmethod = 'nearest';
            source_GA_int              = ft_sourceinterpolate(cfg_source, grandave, template_mri);
            save([save_filename '_interpo'], 'source_GA_int', '-v7.3');
            %}

            % export averaged blob to NifTi format, then use xjview to read out 
            % what brain regions the averaged blob contains
            % (set "intensity" threshold to 2 / 2.5 / 3 (for z-scores), 
            % then press "report" & see console output)
            cfg           = [];
            cfg.filetype  = 'nifti';
            cfg.filename  = save_filename;
            cfg.parameter = 'pow'; %'zscores'
            ft_sourcewrite(cfg, grandave);
        end

        % Alternative ways to read out the anatomical label of a blob: 
        % (1) by looking up an atlas during ft_sourceplot:
        % http://www.fieldtriptoolbox.org/tutorial/aarhus/beamformingerf
        % Doesn't work: (see error below - can't get the Nifti_SPM <-> MNI conversion working)
        %Error using ft_sourceplot:
        %coordinate systems do not match (template mri in Nifti_SPM coords, atlas in MNI coords)
        %
        % (2) use ft_volumelookup:
        % http://www.fieldtriptoolbox.org/faq/how_can_i_determine_the_anatomical_label_of_a_source
    end
    
%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUBFUNCTIONS for source localisation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % calc the difference btwn 2 conditions at all source vertices,
    % so that we can see where in the brain the difference is largest;
    % can also calc z-scores instead of actual power (i.e. which brain region 
    % shows an effect which is much larger than the rest of the brain).
    %
    % @param cond1, cond2: source result for each cond, as obtained from ft_sourceanalysis
    % @param template_mri: for interpolation into common space
    % @param save_filename: where to save 'source_int' (the blob in common space)
    % 
    function localise_effect_source(cond1, cond2, template_sourcemodel, template_mri, save_filename)
        % Compute Percentage Power Difference btwn the 2 conditions
        sourceDiff         = cond2;
        sourceDiff.avg.pow = (cond2.avg.pow - cond1.avg.pow) ./ cond1.avg.pow;

        % replace subject grid with MNI grid (no need to warp, 
        % because the subject-space grid points already have a direct
        % 1-to-1 mapping onto the template grid points)
        sourceDiff.pos = template_sourcemodel.pos;

        % calc z-scores for all vertices
        %{
        x = sourceDiff.avg.pow; % make a shorthand
        sourceDiff.zscores = (x - nanmean(x)) / nanstd(x);
        save([save_filename '_zscores'], 'sourceDiff');
        %}        
        
        % Smear the point solution into an MRI overlay
        % This step does not do any warping; it just adds an mri overlay onto the plot & increases the grid resolution to match template mri (1mm)
        % Q: should we interpolate the blob directly onto template mri, or interpolate onto individual mri first then volumenormalise?
        % A: directly onto template mri, because the subject-space grid points already have a direct 1-to-1 mapping onto the template grid points.
        cfg_source              = [];
        cfg_source.voxelcoord   = 'no';
        cfg_source.parameter    = 'pow';
        cfg_source.interpmethod = 'nearest';
        source_int              = ft_sourceinterpolate(cfg_source, sourceDiff, template_mri);

        save(save_filename, 'source_int'); 

        
        % Plot the result using ft_sourceplot
        %{
        cfg_source              = [];
        %cfg_source.atlas        = fullfile(templates_dir, 'ROI_MNI_V4.nii'); % for looking up the anatomical label
                                                                      % of the source identified
                                                                      % coordsys = 'mni'
        cfg_source.method       = 'ortho';
        cfg_source.funparameter = 'zscores';
        %cfg_source.funcolorlim  = [-0.1 0.1]; % Do this programmatically
        cfg_source.opacitylim   = 'zeromax';
        %cfg_source.location     = [64 -32 8];
        cfg_source.funcolormap  = 'jet';
        ft_sourceplot(cfg_source, source_int);
        %}
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUBFUNCTIONS for ROI analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % retrieve the xyz-coordinates of the given vertices
    function coordinates = get_coordinates_for_vertices(sourcemodel, vertices)
        % each cycle processes one vertex
        for i = 1:length(vertices)
            coordinates(i,:) = sourcemodel.pos(vertices(i), :); %take the warped coordinates
        end

        % alternative method:
        %coordinates = sourcemodel.pos(vertices,:);
    end

    % Load coordinates of centroids from file, then warp into subject space
    % (in preparation for calling create_virtual_sensor_Centroid)
    %
    % @param mri_realigned: correctly realigned MRI file
    % 
    % @output pos_grid: centres of mass in individual space
    %
    function pos_grid = load_centroids(templates_dir, mri_realigned)
        % Load centre of mass information (in mm)
        % These coordinates were generated using:
        % https://github.com/mingruixia/BrainNet-Viewer/blob/master/BrainNet_GenCoord.m
        centre_of_mass = load([templates_dir 'Node_AAL116.txt']);

        % convert mri to mm for consistency
        mri_realigned = ft_convert_units(mri_realigned, 'mm');

        % Align centres of mass into subject space 
        % I.e. we are performing inverse non-linear warping from MNI-->individual
        %cfg = [];
        %cfg.template = fullfile(templates_dir, 'single_subj_T1.nii'); % template brain in MNI space (matches AAL atlas)
        %cfg.nonlinear   = 'yes';
        norm = ft_volumenormalise([], mri_realigned); % transformation matrix from MNI <--> individual
        posback = ft_warp_apply(norm.params, centre_of_mass, 'sn2individual');
        pos_grid = ft_warp_apply(pinv(norm.initial), posback); % xyz-coordinates of each parcel's centre of mass, in individual space

        % convert back to cm
        pos_grid = pos_grid./10; 
    end

    % plot vertices & centroid in the given ROI
    function plot_ROI_centre_of_mass(ROI_name, vertices_coords, centroid, headmodel_singleshell)

        figure('Name',ROI_name, 'NumberTitle','off'); hold on;
        ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor', 'none');
        alpha 0.5; camlight; hold on;

        try           
            ft_plot_mesh(vertices_coords, 'vertexcolor','blue', 'vertexsize',10); hold on;
            ft_plot_mesh(centroid, 'vertexcolor','red', 'vertexsize',30); hold off;
        catch
            fprintf('In ROI "%s": Cannot plot vertices & centroid.\n', ROI_name);
        end
    end

    % create a virtual sensor for the given ROI (based on the spatial filter
    % for each vertex within the ROI), to estimate its activity
    %
    % uses Centroid method for collapsing all vertices within the ROI into a single VE
    %
    % @param vertices:      a list of vertices within this ROI (specified by their vertex index)
    % @param vertices_filters: a list of vertices within this ROI, each entry contains the spatial filter (set of weights) for one vertex
    % @param erf_combined:  average erf across all 4 conditions
    % @param erf:           separate erf for each condition
    %
    function VE = create_virtual_sensor_Centroid(ROI_name, vertices, vertices_filters, erf_combined, erf, conds, headmodel, sourcemodel)

        % compute the coordinates of the centroid
        vertices_coords = get_coordinates_for_vertices(sourcemodel, vertices); % coordinates are in cm
        centroid = find_centroid(vertices_coords); % compute the centroid of this ROI 
        %plot_ROI_centre_of_mass(ROI_name, vertices_coords, centroid, headmodel); % for quality check


        % Create VE by collapsing the activities at all vertices within ROI into one timecourse,
        % weighting the activity at each vertex based on its distance from centre of mass
        %
        %"To generate a single regional timecourse, individual voxel signals
        % were weighted according to their distance from the centre of mass"
        % Brookes et al., (2016)

        VE = [];

        % if there are any vertices in this ROI, create virtual sensor to represent this ROI
        if ~isempty(vertices_filters)
            for i = conds
                % put it into a timelock structure for later calling ft_timelockstatistics (in stats_ROI.m)
                VE.(eventnames_real{i}).time = erf.(eventnames_real{i}).time;
                VE.(eventnames_real{i}).label = {ROI_name};
                VE.(eventnames_real{i}).dimord = 'chan_time';

                % a list of reconstructed source activities, one for each vertex
                timecourses = vertices_filters(:,:) * erf.(eventnames_real{i}).avg(:,:); % estimated source activity = filter * erf (i.e. s = w * X) 

                % weight the timecourse for each vertex by its distance to centre
                timecourses_weighted = [];
                for vertex = 1:size(vertices_coords,1) % each cycle processes one vertex
                    distance = pdist([vertices_coords(vertex,:).*10 ; centroid.*10]); % calculate distance to centre in mm
                    timecourses_weighted(vertex,:) = exp((-distance.^2)./400) .* timecourses(vertex,:); % weight
                end

                % take the mean of all the weighted timecourses (i.e. collapse into one timecourse)
                VE.(eventnames_real{i}).avg = mean(timecourses_weighted, 1);
            end

            % Preserve .sampleinfo field to avoid warnings later
            %VE.sampleinfo = data.sampleinfo;
        end    
    end    

    % create a virtual sensor for the given ROI (based on the spatial filter
    % for each vertex within the ROI), to estimate its activity
    %
    % uses SVD method for collapsing all vertices within the ROI into a single VE
    %
    % @param vertices_filters: a list of vertices within this ROI, each entry contains the spatial filter (set of weights) for one vertex
    % @param erf_combined:  average erf across all 4 conditions
    % @param erf:           separate erf for each condition
    %
    function VE = create_virtual_sensor_SVD(ROI_name, vertices_filters, erf_combined, erf, conds)
        VE = [];

        % if there are any vertices in this ROI, create virtual sensor to represent this ROI
        if ~isempty(vertices_filters)
            % we want to combine the spatial filters for all vertices in this ROI into one filter
            % We do this by performing PCA on concatenated filters * sensor-level cov matrix
            F       = vertices_filters; % make a shorthand
            [u,s,v] = svd(F * erf_combined.cov * F'); % Inner dimensions = number of channels i.e. 125 for KIT child system
            filter  = u' * F;

            % Create VE using this filter
            %VE.label = {ROI_name};
            for i = conds
                % put it into a timelock structure for later calling ft_timelockstatistics (in stats_ROI.m)
                VE.(eventnames_real{i}).time = erf.(eventnames_real{i}).time;
                VE.(eventnames_real{i}).avg(1,:) = filter(1,:) * erf.(eventnames_real{i}).avg(:,:); % estimated source activity = filter * erf (i.e. s = w * X)
                VE.(eventnames_real{i}).label = {ROI_name};
                VE.(eventnames_real{i}).dimord = 'chan_time';
            end

            % Preserve .sampleinfo field to avoid warnings later
            %VE.sampleinfo = data.sampleinfo;


            % only for time-freq analysis
            %{
            %    tfve_plot = figure;
            %    ve_plot   = figure;

            % Create TFR of the VE
            cfg         = [];
            cfg.channel = 'all';
            cfg.method  = 'wavelet';
            cfg.width   = 7;
            cfg.output  = 'pow';
            cfg.foi     = freqband(1):0.25:freqband(2);
            cfg.toi     = -1.0:0.05:1.0; % Need to clean these time windows up for consistency
            cfg.pad     = 'nextpow2';
            TFRwave     = ft_freqanalysis(cfg, eval(ROIs_label{k}));
            save([fname,'_',ROIs_label{k},'_TFRwave_beta'],'TFRwave')

            % Plot
            figure(tfve_plot)
            subplot(4,2,k)
            cfg          = [];
            cfg.ylim     = freqband;
            cfg.baseline = [-0.6 -0.1]; % To do : fix time windows
            cfg.xlim     = [-0.2 1.0];
            ft_singleplotTFR(cfg, TFRwave);
            title(sprintf('%s',ROIs_label{k}));
            colormap(colours)

            figure(ve_plot);
            cfg      = [];
            cfg.xlim = [-0.2 1];
            eval([ROIs_label{k},'_avg_VE_beta           = ft_timelockanalysis(cfg,',ROIs_label{k},');']);
            subplot(4,2,k)
            eval(['ft_singleplotER(cfg,',ROIs_label{k},'_avg_VE_beta)']) % Could add smooth here
            save([fname,'_',ROIs_label{k},'_avg_VE_beta'],[ROIs_label{k},'_avg_VE_beta'])
            %}
        end
    end    

    % create a virtual sensor for the given ROI
    % use this function if you chose cfg.fixedori='no' in ft_sourceanalysis
    % https://mne.tools/stable/auto_tutorials/source-modeling/plot_dipole_orientations.html
    %
    % step (1): compute the reconstructed source activity for each vertex 
    % along all 3 dipole directions (i.e. obtaining 3 timecourses for this vertex)
    % step (2): do a vector combination at every time sample, in order to 
    % combine into 1 timecourse (only keep the vector magnitude, discard the orientation)
    % step (3): collapse all vertices into one VE by taking a plain average
    % over all vertices (the magnitude of activity is positive at all vertices)
    %
    % @param vertices_filters: 3 x channel x vertex (i.e. for each vertex
    %                          in this ROI, 3 sets of weights are provided)
    % @param erf_combined:     average erf across all 4 conditions
    % @param erf:              separate erf for each condition
    %
    function VE = create_virtual_sensor_freeori(ROI_name, vertices_filters, erf_combined, erf, conds)
        VE = [];
                
        % each cycle handles one cond
        for i = conds
            all_vertices_timecourses = [];
            
            % each cycle handles one vertex
            for v = 1:size(vertices_filters, 3)
                % grab the 3 sets of weights for this vertex
                F = vertices_filters(:,:,v);
            
                % compute source timecourse in all 3 dipole orientations
                timecourses_xyz = F * erf.(eventnames_real{i}).avg(:,:); % estimated source activity = filter * erf (i.e. s = w * X)
                timecourse_combined = vecnorm(timecourses_xyz); % vector combination at every time sample, 
                                                             % to obtain the length of the vector (i.e. absolute magnitude 
                                                             % of brain activity, regardless of orientation)
                % add the timecourse for this vertex to the list
                all_vertices_timecourses = [all_vertices_timecourses; timecourse_combined];
            end
            
            % take a plain average over all vertices
            % can also try PCA
            VE_timecourse = mean(all_vertices_timecourses);

            % put it into a timelock structure for later calling ft_timelockstatistics (in stats_ROI.m)
            VE.(eventnames_real{i}).time = erf.(eventnames_real{i}).time;
            VE.(eventnames_real{i}).avg = VE_timecourse;
            VE.(eventnames_real{i}).label = {ROI_name};
            VE.(eventnames_real{i}).dimord = 'chan_time';
        end        
    end    

end % main function end