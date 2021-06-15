% Run ICA to identify eye blinks & other large muscle artefacts 
% (e.g. jaw clenching, hand movements)
%
% @param filter_first - run ICA on 1Hz-filtered data (usually applicable to EEG data,
%                       which has fewer channels, e.g. 32-channel. If you don't filter at 1Hz first,
%                       slow waves will take up too many comps & you might not get a comp for eye blinks)
% if yes: then need to supply these arguments: rawfile, arft, selChLabel
%       e.g. [comp] = ICA_run(true, rawfile, arft, selChLabel);
% if no:  then need to supply this argument: alldata
%       e.g. [comp] = ICA_run(false, alldata);

function [comp] = ICA_run(filter_first, varargin)

    %% prep the data for ICA (apply 1Hz HPF if need to) 
    if (filter_first)
        fprintf(['\nWill run ICA on 1Hz-filtered data\n']); 
        
        rawfile    = varargin{1};
        arft       = varargin{2};
        selChLabel = varargin{3};
        
        hdr = ft_read_header(rawfile, 'dataformat','yokogawa_con'); % read header file
    
        % 1Hz high-pass filter
        cfg                         = [];
        cfg.trialfun                = 'ft_trialfun_general';  
        cfg.channel                 = hdr.label; %hdr.grad.label; (for MEG)
        cfg.continuous              = 'yes';
        cfg.hpfilter                = 'yes';
        cfg.hpfilttype              = 'firws';
        cfg.hpfreq                  = 0.1; % changed to 0.1Hz coz 1Hz actually removed most eye artefacts PRIOR TO ICA (0.5Hz removed quite a lot too)
        cfg.hpfiltdf                = 0.2; % if ICA quality is not good (ie. too many slow drift components), try 0.2-0.3Hz!
        cfg.hpfiltwintype           = 'blackman';
        cfg.hpfiltdir               = 'onepass-zerophase';
        %cfg.dftfreq                 = 50; % removal line noise
        cfg.headerfile              = rawfile;
        cfg.datafile                = rawfile;
        data4ICA                    = ft_preprocessing(cfg);
        
        % low-pass filter (this should be same settings as in filtering.m)
        % if you dont do LPF & bandstop here, you will get ICA comps
        % representing muscle activity & line noise, which cause confusion
        cfg         = [];
        cfg.lpfilter   = 'yes';
        cfg.lpfilttype = 'firws';
        cfg.lpfreq     = 35; % 35 +- 5Hz
        cfg.lpfiltdf   = 10; % wider transition window means it will run much faster
        cfg.lpfiltwintype = 'blackman';
        cfg.lpfiltdir  = 'onepass-zerophase';
        data4ICA = ft_preprocessing(cfg, data4ICA);

        % deal with 50Hz line noise (necessary even after bandpass filter, coz the 50Hz noise is huge)
        cfg          = [];
        cfg.bsfilter = 'yes';
        cfg.bsfreq   = [49.5 50.5];
        % alternatively, can use: (but you need to pad the data to 5~10 seconds)
        % http://www.fieldtriptoolbox.org/faq/what_kind_of_filters_can_i_apply_to_my_data/
        %cfg.dftfilter = 'yes';
        %cfg.dftfreq   = [50 100 150];
        data4ICA = ft_preprocessing(cfg, data4ICA);
        

        % reject the artifacts and channels that were marked in Steps 3 & 4 
        arft.artfctdef.reject       = 'partial';
        data4ICA                    = ft_rejectartifact(arft, data4ICA);

        cfg                         = [];
        cfg.channel                 = selChLabel;
        data4ICA                    = ft_selectdata(cfg, data4ICA);

    else % if not applying the 1Hz filter, just read in the data
        
        data4ICA = varargin{1};
    
    end        

    
    %% Run ICA (don't use 'fastica' method, use 'runica' method)
    disp('About to run ICA using the Runica method')
    cfg            = [];
    cfg.method     = 'runica';
    cfg.channel    = 'all'; 
    comp           = ft_componentanalysis(cfg, data4ICA);

end