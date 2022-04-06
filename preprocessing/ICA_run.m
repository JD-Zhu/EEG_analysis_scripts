% Run ICA to identify eye blinks & other large muscle artefacts 
% (e.g. jaw clenching, hand movements)
%
% @param filter_first - do another HPF before running ICA (usually applicable to EEG data,
%                       which has fewer channels, e.g. 32-channel. If you don't filter at 1Hz first,
%                       slow waves could take up too many comps & you might not get a good comp for eye blinks)
% if yes: then need to supply these arguments: filters (hpfreq, hpfiltdf, lpfreq, lpfiltdf), rawfile, arft, selChLabel
%       e.g. [comp] = ICA_run(true, [0.1 0.2 35 10], rawfile, arft, selChLabel);
% if no:  then need to supply this argument: alldata
%       e.g. [comp] = ICA_run(false, alldata);

function [comp] = ICA_run(filter_first, varargin)

    %% prep the data for ICA (apply another HPF if need to) 
    if (filter_first)
        filters         = varargin{1};
        rawfile         = varargin{2};
        arft            = varargin{3};
        selChLabel      = varargin{4};
        
        fprintf(['\nWill run ICA on ' num2str(filters(1)) 'Hz-filtered data\n']); 
        
        hdr = ft_read_header(rawfile, 'dataformat','yokogawa_con'); % read header file
    
        % 1Hz high-pass filter
        cfg                         = [];
        cfg.trialfun                = 'ft_trialfun_general';  
        cfg.channel                 = hdr.label; %hdr.grad.label; (for MEG)
        cfg.continuous              = 'yes';
        cfg.hpfilter                = 'yes';
        cfg.hpfilttype              = 'firws';
        cfg.hpfreq                  = filters(1);
        cfg.hpfiltdf                = filters(2);
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
        cfg.lpfreq     = filters(3);
        cfg.lpfiltdf   = filters(4);
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
        

        % For NeuroPrax EEG, remove the prefix "EEG " from channel labels
        % (to be consistent with alldata - we did this in main.m)
        if contains(data4ICA.label(1), 'EEG')
            temp = data4ICA.label(1:27);
            temp = cellfun(@(x) x(5:end), temp, 'un', 0); % remove first 4 chars in each cell
            data4ICA.label(1:27) = temp;
            data4ICA.hdr.label(1:27) = temp; % also update this (just in case)
        end

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