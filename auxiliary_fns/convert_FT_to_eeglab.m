% Convert data from fieldtrip format to eeglab format (i.e. raw matrix), 
% so that they can be accepted by various stats toolbox (e.g. TFCE, EMSf)
%
%    'eeglab' format = 3d matrix of values (e.g. subj x chan x time)
% 'fieldtrip' format = cell array of subjects/trials
%
% @param data:          can be ERF result or ROI result
% @param which_toolbox: 'TFCE' or 'EMSf' (diff toolboxes rq diff formats)
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function subj_chan_time = convert_FT_to_eeglab(data, which_toolbox)
    subj_chan_time = []; % initialise the 3d matrix

    for subject = 1:length(data) % loop thru all subjects/trials
        % add this cell's "chan x time" matrix to the 3d matrix
        if strcmp(which_toolbox, 'TFCE')
            chan_time = data{subject}.avg; % each cell contains the timelock struct for one subject (.avg field contains chan x time matrix)
        elseif strcmp(which_toolbox, 'EMSf')
            chan_time = data{subject}; % each cell contains the chan x time matrix for one trial
        else
            error('Error in convert_FT_to_eeglab(): the toolbox specified is invalid.');
        end
        subj_chan_time = cat(3, subj_chan_time, chan_time); % concatenate along the 3rd dimension
    end
    
    % subj_chan_time is now the 3d matrix containing all subjects ("subject" being the 3rd dimension)
    % if using 'TFCE' toolbox, change the order of matrix dimensions to: subj x chan x time
    if strcmp(which_toolbox, 'TFCE')
        subj_chan_time = permute(subj_chan_time, [3 1 2]); 
    end
end