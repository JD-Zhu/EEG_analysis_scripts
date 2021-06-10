% wrapper function for calling ept_TFCE(), so that settings only need to be changed in one place
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function Results = myWrapper_ept_TFCE(data1, data2)

    % the input format to TFCE requires a "channel" dimension, so we check this
    % https://github.com/Mensen/ept_TFCE-matlab/issues/21
    if (size(data1, 2) == 1) % if there is only 1 channel (e.g. ROI time series),
                             % we fake a 2nd channel by making a copy of the 1st channel
        data1 = repmat(data1, [1,2,1]);
        data2 = repmat(data2, [1,2,1]);

        % we also need to treat this data as time-frequency data (so that TFCE won't look for a channel locations file)
        flag_ft = true;
        chanlocs = [];

    else % for normal (multi-channel) data, e.g. ERF time series
        
        flag_ft = false;
        
        % read the channel locations
        %addpath(genpath('H:\eeglab_current\eeglab14_1_1b\'));
        %chanlocs = readlocs('chanlocs_XYZ.txt', 'filetype','custom', 'format',{'X','Y','Z'});
        chanlocs = []; load('chanlocs.mat');

    end       

    
    Results = ept_TFCE(data1, data2, ...
        chanlocs, ...
        'type', 'd', ...
        'flag_ft', flag_ft, ...
        'flag_tfce', true, ...
        'nPerm', 2000, ...
        'rSample', 200, ...
        'flag_save', false);
        %'saveName', [ResultsFolder 'TFCE_temp\\ept_' ROI_name '.mat']); % set a location to temporarily store the output. we don't need to save it, but if you don't set a location, it will litter arond your current directory

end