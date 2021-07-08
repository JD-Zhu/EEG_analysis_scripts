% do_HPF: if yes, apply 0.01Hz high-pass filter

function [alldata] = preprocessing(rawfile, do_HPF)
    
    hdr = ft_read_header(rawfile);%, 'dataformat','yokogawa_con'); % read header file
         
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

    %{
    cfg.bpfilter   = 'yes';
    cfg.bpfreq     = [0.2 30]; % bandpass filter [0.5 30], successfully filtered out the low-freq drift!!
    cfg.bpfiltwintype = 'blackman'; % may help to get rid of the ringing effect?
    cfg.bpfiltord = 64; % default is 4
    alldata = ft_preprocessing(cfg);
    %}        

    % highpass filter (optional)
    if (do_HPF)
        cfg            = [];
        cfg.hpfilter   = 'yes';
        cfg.hpfilttype = 'firws';
        cfg.hpfreq     = 0.01; % 0.01 +- 0.01Hz
        cfg.hpfiltdf   = 0.02; % transition window width (for firws; this param overrides order)
                              % hpfreq - (hpfiltdf / 2) must be >= 0
        cfg.hpfiltwintype = 'blackman';
        cfg.hpfiltdir  = 'onepass-zerophase';
        alldata = ft_preprocessing(cfg, alldata);
    end

    % lowpass filter
    cfg         = [];
    cfg.lpfilter   = 'yes';
    cfg.lpfilttype = 'firws';
    cfg.lpfreq     = 35; % 35 +- 5Hz
    cfg.lpfiltdf   = 10; % wider transition window means it will run much faster
    cfg.lpfiltwintype = 'blackman';
    cfg.lpfiltdir  = 'onepass-zerophase';
    alldata = ft_preprocessing(cfg, alldata);

    % deal with 50Hz line noise (necessary even after bandpass filter, coz the 50Hz noise is huge)
    cfg          = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = [49.5 50.5];
    % alternatively, can use: (but you need to pad the data to 5~10 seconds)
    % http://www.fieldtriptoolbox.org/faq/what_kind_of_filters_can_i_apply_to_my_data/
    %cfg.dftfilter = 'yes';
    %cfg.dftfreq   = [50 100 150];
    alldata = ft_preprocessing(cfg, alldata);
    
end