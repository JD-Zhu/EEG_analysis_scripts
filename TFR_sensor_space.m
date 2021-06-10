function TFR_sensor_space(rawfile,do_plot)

hdr = ft_read_header(rawfile,'dataformat','yokogawa_con'); %hdr to get Fs etc.

if ~isempty(find(ismember(hdr.label,'AG160')))
    fprintf('ADULT MEG\n')
    isadult = 1;
else
    isadult = 0;
end
% trigger events
eventcodes = {{'cue'},{'1'};{'imperative'},{'2'};{'response'},{'3'}};

% ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
cfg                      = [];
cfg.headerfile           = rawfile;
cfg.datafile             = rawfile;
cfg.trialdef.triallength = Inf;
cfg.trialdef.ntrials     = 1; % read in all data as a single segment
cfg                      = ft_definetrial(cfg);

% ft_preprocessing: reads in MEG data
cfg.continuous = 'yes';
cfg.bpfilter   = 'yes';
cfg.bpfreq     = [1.0 40]; % bandpass filter

alldata = ft_preprocessing(cfg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%here we can convert to child MEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isadult
    % Select channels 1-125 (i.e. MEG data)
    cfg           = [];
    cfg.channel   = alldata.label(1:160);
    alldata       = ft_selectdata(cfg, alldata);
    [alldata,~,~] = convert_adult_2_child(alldata);
else
    % Select channels 1-125 (i.e. MEG data)
    cfg         = [];
    cfg.channel = alldata.label(1:125);
    alldata     = ft_selectdata(cfg, alldata);
end

% deal with 50Hz line noise
cfg          = [];
cfg.bsfilter = 'yes';
cfg.bsfreq   = [49.5 50.5];
alldata      = ft_preprocessing(cfg, alldata);

%     Create layout file for later + save
cfg      = [];
cfg.grad = alldata.grad; % struct containing gradiometer definition
lay      = ft_prepare_layout(cfg, alldata); % creates a 2-D layout of the channel locations

load([cd,'/event.mat']);
load([cd,'/trial.mat']);

trialinfo_b = event;

cfg     = [];
cfg.trl = trl;
alldata = ft_redefinetrial(cfg, alldata);

cfg         = [];
cfg.demean  = 'yes';
cfg.detrend = 'yes';

alldata = ft_preprocessing(cfg, alldata);

eventnames = eventcodes(:,1); % extract a list of all event names
eventnames = [eventnames{:}]; % convert into strings
trialsgone = 0;

% each cycle is one "block" (i.e. one '.con' file)
for i = 1:length(trialinfo_b)
    for j = 1:length(eventcodes)
        events.(eventnames{j}) = find(strcmp({trialinfo_b(i).value}, eventcodes{j,2})); % 9 fields representing the 9 types of events
        % each field contains a list of all events belonging to this type
        % (by matching event code)
        % NB. this is a tmp var, it gets overwritten in each cycle
    end
    if i == 1 % first block only
        for j = 1:length(eventcodes)
            events_b.(eventnames{j}) = events.(eventnames{j})'; % save the lists to a perm var, also transpose each list
        end
    else % all other blocks
        trialsinblock = length(trialinfo_b(i-1)); % how many "trials" (i.e. events) were identified in previous block
        trialsgone    = trialsgone + trialsinblock; % add this number to the total number of "past" trials
        
        for j = 1:length(eventcodes)
            events_b.(eventnames{j}) = [events_b.(eventnames{j}); events.(eventnames{j})' + trialsgone]; % continue to append to the perm lists (stored in "events_b")
            % in the end, each perm list will contain all events of that type from all blocks
        end
    end
end

% ft_redefine the "response" trials & calc erf
cfg        = [];
cfg.trials = events_b.response; % list of "response" events
cfg.toilim = [-0.5 0.5];
response   = ft_redefinetrial(cfg, alldata);

cfg          = [];
response_erf = ft_timelockanalysis(cfg, response);

%Run ICA on the "response" trials
disp('About to run ICA using the SVD method')
cfg           = [];
cfg.method    = 'svd';
response_comp = ft_componentanalysis(cfg, response_erf);

%Change the colourmap
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colours = colormap(flipud(brewermap(64, 'RdBu'))); % change the colormap

if do_plot==1
    %Display Components - change layout as needed
    cfg                 = [];
    cfg.viewmode        = 'component';
    cfg.layout          = lay;
    cfg.channel         = {'svd001','svd002','svd003'};
    cfg.channelcolormap = colours;
    ft_databrowser(cfg, response_comp)
else
end

%**unmix here...is that necessary?***
cfg           = [];
cfg.component = 1:3; % reject top 5 components
alldata_clean = ft_rejectcomponent(cfg, response_comp, alldata); % reject these comps from all trials

% Reject Outlier Trials
[z,bad_trials,alldata_clean] = artifacts_max_z(alldata_clean,10);

% ft_redefine all other event types & calc erf
for j = 1:length(eventcodes)
    cfg                    = [];
    cfg.trials             = events_b.(eventnames{j});
    trials.(eventnames{j}) = ft_redefinetrial(cfg, alldata);
end

% do the same for the cleaned data
for j = 1:length(eventcodes)
    cfg                          = [];
    cfg.trials                   = events_b.(eventnames{j});
    trials_clean.(eventnames{j}) = ft_redefinetrial(cfg, alldata_clean);
end

% calc erf
for j = 1:length(eventcodes)
    cfg                 = [];
    cfg.nanmean         = 'yes';
    erf.(eventnames{j}) = ft_selectdata(cfg, trials.(eventnames{j})); % Do this because we kept bad trials as NaN
    erf.(eventnames{j}) = ft_timelockanalysis(cfg, erf.(eventnames{j})); % Do this to create average field and timelock struct
end

% fill in some erf struct
for j = 1:length(eventcodes)
    tmp                        = trials.(eventnames{j}).trial;
    tmp_mean                   = nanmean(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),size(tmp,2)),3);
    erf.(eventnames{j}).avg    = tmp_mean;
    erf.(eventnames{j}).time   = erf.(eventnames{j}).time;
    erf.(eventnames{j}).dimord = 'chan_time';
end

% do the same for the cleaned data
for j = 1:length(eventcodes)
    tmp                              = trials_clean.(eventnames{j}).trial;
    tmp_mean                         = nanmean(reshape([tmp{:}],size(tmp{1},1),size(tmp{1},2),size(tmp,2)),3);
    erf_clean.(eventnames{j})        = trials_clean.(eventnames{j});
    erf_clean.(eventnames{j}).avg    = tmp_mean;
    erf_clean.(eventnames{j}).time   = erf_clean.(eventnames{j}).time{1};
    erf_clean.(eventnames{j})        = rmfield(erf_clean.(eventnames{j}), 'trial');
    erf_clean.(eventnames{j}).dimord = 'chan_time';
end

cfg              = [];
cfg.showlabels   = 'yes';
cfg.fontsize     = 6;
cfg.layout       = lay;
cfg.baseline     = [-0.2 0];
cfg.baselinetype = 'absolute';

if do_plot==1
    for j = 1:length(eventcodes)
        figure;
        ft_multiplotER(cfg, erf.(eventnames{j}), erf_clean.(eventnames{j}));
    end
else
end

% to compare the 3 conds
% (requires the list of 'cfg' assignments above, if running in console)
figure('Name','ft_multiplotER: erf_clean.cue, erf_clean.imperative, erf_clean.response');
ft_multiplotER(cfg, erf_clean.cue, erf_clean.imperative, erf_clean.response);

% Calc global averages across all sensors (GFP = global field potentials)
cfg        = [];
cfg.method = 'power';
for j = 1:length(eventcodes)
    erf_clean_GFP.(eventnames{j}) = ft_globalmeanfield(cfg, erf_clean.(eventnames{j}));
end

% plot global averages for cue-locked
cfg = [];
figure('Name','cue_GFP'); hold on
for j = 1:3
    plot(erf_clean_GFP.(eventnames{j}).time, erf_clean_GFP.(eventnames{j}).avg);
end
legend([eventcodes{1:3, 1}])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TF analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = trials_clean.cue;

alpha = [8 15];
beta  = [16 30];
mu    = [7.5 12.5];

cfg               = [];
cfg.method        = 'distance';
cfg.feedback      = 'no'; % show a neighbour plot
cfg.neighbourdist = 5;
cfg.grad          = data.grad;
neighbours        = ft_prepare_neighbours(cfg);

cfg              = [];
cfg.feedback     = 'no';
cfg.method       = 'template';
cfg.planarmethod = 'sincos';
cfg.channel      = {'MEG'};
cfg.trials       = 'all';
cfg.neighbours   = neighbours;
data_planar      = ft_megplanar(cfg,data);

%% time-frequency analysis on axial and synthetic planar data
cfg            = [];
cfg.output     = 'pow';
cfg.channel    = 'MEG';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.foi        = 2:1:40;
cfg.toi        = [-1:0.05:1];
cfg.trials     = 'all';
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
cfg.pad        = 'nextpow2';
cfg.keeptrials = 'yes';%% we keep the single trials
freq_planar    = ft_freqanalysis(cfg,data_planar);

cfg               = [];
cfg.combinemethod = 'sum';
freq_combined     = ft_combineplanar(cfg,freq_planar);
TFRwave           = ft_freqdescriptives([],freq_combined);

cfg           = [];
cfg.latency   = [-0.3 inf];
cfg.frequency = alpha;
TFRwave_alpha = ft_selectdata(cfg,TFRwave);
cfg.frequency = beta;
TFRwave_beta  = ft_selectdata(cfg,TFRwave);
cfg.frequency = mu;
TFRwave_mu    = ft_selectdata(cfg,TFRwave);

cfg              = [];
cfg.baseline     = [-0.3 -0.1];
cfg.baselinetype = 'relchange';

norm_TFR_broad = ft_freqbaseline(cfg,TFRwave);
norm_TFR_alpha = ft_freqbaseline(cfg,TFRwave_alpha);
norm_TFR_beta  = ft_freqbaseline(cfg,TFRwave_beta);
norm_TFR_mu    = ft_freqbaseline(cfg,TFRwave_mu);

patches1 = {neighbours.neighblabel};
patches2 = {neighbours.label}';
for i=1:length(patches1)
    patches{i} = vertcat(patches1{i},patches2(i)); %add the seed sensor into the patch
end


%%%%
%ERD
%%%%%
minima_alpha                                   = (ismember(norm_TFR_alpha.powspctrm,min(min(min(norm_TFR_alpha.powspctrm(:,:,:))))));
[chan_alpha_erd freq_alpha_erd time_alpha_erd] = ind2sub(size(minima_alpha),find(minima_alpha));

minima_beta                                 = (ismember(norm_TFR_beta.powspctrm,min(min(min(norm_TFR_beta.powspctrm(:,:,:))))));
[chan_beta_erd freq_beta_erd time_beta_erd] = ind2sub(size(minima_beta),find(minima_beta));

minima_mu                             = (ismember(norm_TFR_mu.powspctrm,min(min(min(norm_TFR_mu.powspctrm(:,:,:))))));
[chan_mu_erd freq_mu_erd time_mu_erd] = ind2sub(size(minima_mu),find(minima_mu));

%%%%%%
%ERS
%%%%%%
maxima_alpha                                   = (ismember(norm_TFR_alpha.powspctrm,max(max(max(norm_TFR_alpha.powspctrm(:,:,:))))));
[chan_alpha_ers freq_alpha_ers time_alpha_ers] = ind2sub(size(maxima_alpha),find(maxima_alpha));

maxima_beta                                 = (ismember(norm_TFR_beta.powspctrm,max(max(max(norm_TFR_beta.powspctrm(:,:,:))))));
[chan_beta_ers freq_beta_ers time_beta_ers] = ind2sub(size(maxima_beta),find(maxima_beta));

maxima_mu                             = (ismember(norm_TFR_mu.powspctrm,max(max(max(norm_TFR_mu.powspctrm(:,:,:))))));
[chan_mu_ers freq_mu_ers time_mu_ers] = ind2sub(size(maxima_mu),find(maxima_mu));

cfg      = [];
cfg.xlim = [-0.3 1];
cfg.ylim = [5 30];

if do_plot==1
    %%%
    %ERD
    %%%
    figure;
    
    subplot(3,2,1)
    cfg.channel  = patches{chan_alpha_erd}; % top figure
    cfg.colormap = colours;
    ft_singleplotTFR(cfg, norm_TFR_broad);
    title('Alpha')
    
    subplot(3,2,3)
    cfg.channel  = patches{chan_beta_erd}; % top figure
    cfg.colormap = colours;
    ft_singleplotTFR(cfg, norm_TFR_broad);
    title('Beta')
    
    subplot(3,2,5)
    cfg.channel  = patches{chan_mu_erd}; % top figure
    cfg.colormap = colours;
    ft_singleplotTFR(cfg, norm_TFR_broad);
    title('Mu')
    
    cfg                 = [];
    cfg.lay             = lay;
    cfg.highlight       = 'on';
    cfg.highlightsymbol = '*';
    cfg.highlightsize   = 10;
    
    subplot(3,2,2)
    cfg.highlightchannel = data.label(find(ismember({neighbours.label},patches{chan_alpha_erd})));
    cfg.xlim             = [norm_TFR_alpha.time(time_alpha_erd) norm_TFR_alpha.time(time_alpha_erd)];
    cfg.colormap         = colours;
    ft_topoplotTFR(cfg, norm_TFR_alpha)
    
    subplot(3,2,4)
    cfg.highlightchannel = data.label(find(ismember({neighbours.label},patches{chan_beta_erd})));
    cfg.xlim             = [norm_TFR_beta.time(time_beta_erd) norm_TFR_alpha.time(time_beta_erd)];
    cfg.colormap         = colours;
    ft_topoplotTFR(cfg, norm_TFR_beta)
    
    subplot(3,2,6)
    cfg.highlightchannel = data.label(find(ismember({neighbours.label},patches{chan_mu_erd})));
    cfg.xlim             = [norm_TFR_mu.time(time_mu_erd) norm_TFR_alpha.time(time_mu_erd)];
    cfg.colormap         = colours;
    ft_topoplotTFR(cfg, norm_TFR_mu)
    
    %%%%%%%%%%%%%%%%
    %ERS
    %%%%%%%%%%%%%%%%
    cfg      = [];
    cfg.xlim = [-0.3 1];
    cfg.ylim = [5 30];
    
    
    subplot(3,2,1)
    cfg.channel  = patches{chan_alpha_ers}; % top figure
    cfg.colormap = colours;
    ft_singleplotTFR(cfg, norm_TFR_broad);
    title('Alpha')
    
    subplot(3,2,3)
    cfg.channel  = patches{chan_beta_ers}; % top figure
    cfg.colormap = colours;
    ft_singleplotTFR(cfg, norm_TFR_broad);
    title('Beta')
    
    subplot(3,2,5)
    cfg.channel  = patches{chan_mu_ers}; % top figure
    cfg.colormap = colours;
    ft_singleplotTFR(cfg, norm_TFR_broad);
    title('Mu')
    
    cfg                 = [];
    cfg.lay             = lay;
    cfg.highlight       = 'on';
    cfg.highlightsymbol = '*';
    cfg.highlightsize   = 10;
    
    subplot(3,2,2)
    cfg.highlightchannel = data.label(find(ismember({neighbours.label},patches{chan_alpha_ers})));
    cfg.xlim             = [norm_TFR_alpha.time(time_alpha_ers) norm_TFR_alpha.time(time_alpha_ers)];
    cfg.colormap         = colours;
    ft_topoplotTFR(cfg, norm_TFR_alpha)
    
    subplot(3,2,4)
    cfg.highlightchannel = data.label(find(ismember({neighbours.label},patches{chan_beta_ers})));
    cfg.xlim             = [norm_TFR_beta.time(time_beta_ers) norm_TFR_alpha.time(time_beta_ers)];
    cfg.colormap         = colours;
    ft_topoplotTFR(cfg, norm_TFR_beta)
    
    subplot(3,2,6)
    cfg.highlightchannel = data.label(find(ismember({neighbours.label},patches{chan_mu_ers})));
    cfg.xlim             = [norm_TFR_mu.time(time_mu_ers) norm_TFR_alpha.time(time_mu_ers)];
    cfg.colormap         = colours;
    ft_topoplotTFR(cfg, norm_TFR_mu)
    
    %%%%%%%%%%%%%%%%
    %to do calc power over average time 200-800 then find min max etc.
    %%%%%%%%%%%%%%%%
    figure;
    cfg          = [];
    cfg.lay      = lay;
    cfg.xlim     = [0.2 0.8];
    cfg.colormap = colours;
    
    subplot(3,1,1)
    ft_topoplotTFR(cfg, norm_TFR_alpha)
    title('Alpha')
    
    subplot(3,1,2)
    ft_topoplotTFR(cfg, norm_TFR_beta)
    title('Beta')
    
    subplot(3,1,3)
    ft_topoplotTFR(cfg, norm_TFR_mu)
    title('Mu')
    
else
end


save ('TFR', 'norm_TFR_broad')
save ('erf', 'erf')
save ('erf_clean', 'erf_clean')
save ('trials', 'trials')
save ('trials_clean', 'trials_clean')
end