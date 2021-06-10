% NOTE: this function is no longer used!
% shouldn't do any detrending on epoched data, as this causes edge effects ("bowtie" shape)
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function [alldata_clean] = detrend_on_epoched_data(alldata, chanlocs)

    % <Opt 1> FT_detrend
    
    cfg = [];
    cfg.detrend = 'yes'; % removes low-frequency drift
    cfg.polyremoval = 'yes'; % == remove higher order polynomial (default: 2nd-order)
    cfg.polyorder = 40; % customise the order (higher order => affects higher freqs?)

    alldata_clean = ft_preprocessing(cfg, alldata);

    
    % <Opt 2> locdetrend

    % first, convert to EEGLAB format
    data = alldata;
    output_file = 'epoched.mat';
    save(output_file, 'data', '-v7.3');
    [EEG] = fieldtrip2eeglab(output_file);
    
    % fill in the chanlocs var
    EEG.chanlocs = chanlocs;

    %EEG = robust_locdetrend(EEG, 'no', 300, 150, 'no'); % ~0.83Hz HPF; works perfectly (epoches are completely flat after lcodetrend), but is this removing too much?
    EEG = robust_locdetrend(EEG, 'no', 600, 300, 'no'); % ~0.42Hz HPF; a little bit of change to the waveform, still some drifts remaining
    %EEG = robust_locdetrend(EEG, 'no', 900, 450, 'yes'); % ~0.28Hz HPF; a little bit of change to the waveform, still some drifts remaining
    %EEG = robust_locdetrend(EEG, 'no', 1000, 500, 'yes'); % ~0.25Hz HPF; almost no change to the waveform

    % convert back to FT format
    data_temp = eeglab2fieldtrip(EEG, 'preprocessing');  
    
    alldata_clean = alldata;
    alldata_clean.trial = data_temp.trial;     
end