% browse raw data to remove noisy segments (e.g. clenched jaw)
% http://www.fieldtriptoolbox.org/walkthrough/#visual-data-inspection

function [arft] = mark_artefact(alldata)

    cfg           = [];
    cfg.viewmode  = 'vertical';
    cfg.continous = 'yes';

    cfg.blocksize = 120; % display 2-min segments
    cfg.ylim      = [ -256   256 ];
    %{
    % if trigger channels were retained in this data, first use 
    % this scaling to mark the break periods (ie. when there are no triggers)
    if (length(alldata.label) > 160)
        cfg.blocksize = 300; % display 5-min segments
        cfg.ylim      = [ -0.25  0.25 ];
        cfg.channel   = alldata.label(end-1:end); % display channels [191 193];
        
        % Note: after marking the break periods, still need to manually 
        % enter the scaling below, to mark artefact in MEG channels!!
        
    elseif (length(alldata.label) == 160)    
        % this is the best scaling to see artefact in MEG channels
        cfg.blocksize = 120; % display 2-min segments
        cfg.ylim      = [ -4e-13   4e-13 ];
        
    else % should never be here
        fprintf('Error in mark_artefact(): number of MEG channels < 160.\n');        
    end
    %}
    
    % randomly sample a few channels, to check if the bandpass filter 
    % can successfully filter out the low-freq drift
    %{
    %cfg.channel = [1 8 15 23 30 37 44 51 58 65 72 78 84 91 98 104 111 119 125 132 139 146 153]; 
    cfg.channel = [1 11 35 49 60 72 85 98 111 125 139 153]; 
    cfg.ylim = [ -1.0941e-12  1.0941e-12 ];
    %}
    
    arft = ft_databrowser(cfg, alldata);
end