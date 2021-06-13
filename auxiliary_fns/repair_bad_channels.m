% repairs bad/missing/rejected channels after ft_rejectvisual
% can also run directly on the computed erf of each subject
% 
% While ft_rejectvisual provides the option to repair channels 
% straight away after rejection (by setting cfg.keepchannel = 'repair'),
% that approach doesn't allow you to select the 'interpolation' method 
% (it uses the default 'weighted' method, which is no good for bad channels 
% that lie next to each other)
%
% This function assumes the rejected channels have been removed from 
% the channel list (data.label). This would be the case 
% if you set cfg.keepchannel = 'no' when calling ft_rejectvisual.
% The removed channels will be automatically detected in this script.
%
% @param data:       the data with missing channels to be repaired (can be indi trials or final erf)
% @param neighbours: the neighbours file created by ft_prepare_neighbours
% @param all_labels: full list of all channels
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function data = repair_bad_channels(data, neighbours, all_labels)
    % reject all potentially problematic channels first,
    % so that the noise doesn't leak into the repaired channels thru interpolation
    %{
    cfg = [];
    cfg.channel = {'all', '-AG083', '-AG087', '-AG088', '-AG082', '-AG084', '-AG086', '-AG081', '-AG085', '-AG089'};
    data = ft_selectdata(cfg, data); % if any of these channels have already been removed (e.g. during visual rejection), then this won't do anything to those channels
    %}
    
    
    % Ref:
    % %http://www.fieldtriptoolbox.org/example/fixing_a_missing_sensor
    
    % automatically extract the list of missing channels
    % (i.e. the bad channels you manually selected during visual rejection)
    [notmissing, ~] = match_str(all_labels, data.label); % find all retained channels
    goodchans   = false(numel(all_labels),1); % initialise logical array with 0s
    goodchans(notmissing) = true; % mark all good channels as 1
    badchanindx = find(goodchans==0); % the rest are bad channels
    badChannels = all_labels(badchanindx); % get the channel names

    cfg = [];
    cfg.method         = 'spline'; % use spline method for missing channels that lie next to each other
    cfg.badchannel     = badChannels;
    cfg.missingchannel = {};
    cfg.neighbours     = neighbours;
    data = ft_channelrepair(cfg, data);    
end
