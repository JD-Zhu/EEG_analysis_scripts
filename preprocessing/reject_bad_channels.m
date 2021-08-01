% manually identify bad channels
% Note: this fn does not change the data (ie. it does not remove the bad channels), 
% it only produces plots to help us determine which channels are bad

function [selChLabel] = reject_bad_channels(alldata)
    % plot all channels to get an overall idea
    cfg          = [];
    cfg.method   = 'trial';
    ft_rejectvisual(cfg, alldata);
    
    % multiplot 
    %{
    cfg = [];
    cfg.layout = lay;
    figure; ft_multiplotER(cfg, alldata);
    %}            
    
    % plot each channel individually, for manual rejection
    % Note - This is absly not informative! (even flat channels show up
    %normal... prob due to diff scaling?)
    % So, instead of using this, just enter the channels you want to reject
    % in the visual summary (below)
    %{
    cfg                         = [];
    cfg.method                  = 'channel';
    %cfg.method                  = 'summary';
    %cfg.alim                    = 1e-10;
    cfg.keepchannel             = 'no';
    cfg.keeptrial               = 'nan';    
    alldata = ft_rejectvisual(cfg, alldata); 
    %}
    
    % Display visual summary to find outlier channels
    cfg              = [];
    cfg.feedback     = 'no'; % suppress console output (so that it's easy to find the SubjectID we printed out above)
    cfg.method       = 'summary';
    cfg.metric       = 'zvalue'; % default is 'var'
    cfg.keepchannel  = 'no';
    cfg.keeptrial    = 'nan'; % we keep the rejected trials as 'NaN' here,
        % because if we remove them, that will change the indices of all subsequent trials,
        % which will no longer match the indices we are using in events_allBlocks
    alldata = ft_rejectvisual(cfg, alldata);

    
    % return a list of channels to keep
    selChLabel = alldata.label;
    
end
