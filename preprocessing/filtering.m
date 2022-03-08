% @param filters - should be supplied in this format: [hpfreq hpfiltdf lpfreq lpfiltdf]
%                  e.g. [0.01 0.02 35 10] means 0.01+-0.01Hz HPF and 35+-5Hz LPF
%                  Note: hpfreq - (hpfiltdf / 2) must be >= 0
% @param do_HPF - if no, ignore the hpf settings

function [alldata] = filtering(alldata, do_HPF, filters)
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
        cfg.hpfreq     = filters(1);
        cfg.hpfiltdf   = filters(2); % transition window width (for firws; this param overrides order)
        cfg.hpfiltwintype = 'blackman';
        cfg.hpfiltdir  = 'onepass-zerophase';
        alldata = ft_preprocessing(cfg, alldata);
    end

    % lowpass filter
    cfg         = [];
    cfg.lpfilter   = 'yes';
    cfg.lpfilttype = 'firws';
    cfg.lpfreq     = filters(3);
    cfg.lpfiltdf   = filters(4); % wider transition window means it will run much faster
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