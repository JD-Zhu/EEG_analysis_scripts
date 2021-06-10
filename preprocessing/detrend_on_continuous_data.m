% perform detrending on continuous data, using locdetrend()
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function [alldata] = detrend_on_continuous_data(alldata, chanlocs)

    % 1. Convert to EEGLAB format

    data = alldata;
    % remove the trigger channels from the data
    cfg         = [];
    cfg.channel = data.label(1:160);
    data     = ft_selectdata(cfg, data);
    
    tmp_save_file = [pwd '\\..\\data4eeglab.mat'];
    save(tmp_save_file, 'data', '-v7.3');

    eeglab; % start eeglab (which will autoly add paths, so we can call the fn below)
    [EEG] = fieldtrip2eeglab(tmp_save_file);

    % fill in the chanlocs var
    EEG.chanlocs = chanlocs;

    
    % 2. Perform local detrend
    
    % continuous data is currently not supported, so we pretend this to be epoched data
    % (by creating 2 identical trials)
    EEG.data = repmat(EEG.data, [1,1,2]);

    %EEG_detrended = robust_locdetrend(EEG, 'no', 600, 300, 'no'); 
    EEG_detrended = robust_locdetrend(EEG, 'no', 1000, 500, 'no'); % both work equally well, so use the wider window (safer)

    % remove the fake trial we created b4
    EEG_detrended.data = EEG_detrended.data(:,:,1);

    
    % 3. Convert back to FT format
    data_temp = eeglab2fieldtrip(EEG_detrended, 'preprocessing');  
    
    alldata.trial = data_temp.trial;  
    alldata.label = alldata.label(1:160); % make channels consistent
    
end