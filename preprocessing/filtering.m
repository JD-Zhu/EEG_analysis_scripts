% do_HPF: if yes, apply 0.01Hz high-pass filter

function [alldata] = preprocessing(alldata, do_HPF, hpfreq, hpfiltdf, lpfreq, lpfiltdf)
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
        cfg.hpfreq     = hpfreq; % 0.01 +- 0.01Hz
        cfg.hpfiltdf   = hpfiltdf; % transition window width (for firws; this param overrides order)
                              % hpfreq - (hpfiltdf / 2) must be >= 0
        cfg.hpfiltwintype = 'blackman';
        cfg.hpfiltdir  = 'onepass-zerophase';
        alldata = ft_preprocessing(cfg, alldata);
    end

    % lowpass filter
    cfg         = [];
    cfg.lpfilter   = 'yes';
    cfg.lpfilttype = 'firws';
    cfg.lpfreq     = lpfreq; % 35 +- 5Hz
    cfg.lpfiltdf   = lpfiltdf; % wider transition window means it will run much faster
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