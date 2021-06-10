% Remove trials containing NaNs (i.e. any trials containing the sections
% we manually marked previously)
function [all_blocks_clean, clean_event] = remove_nan_trials(all_blocks, event)
    % check for boundary case - when there are no trials in the input data
    if isempty(event)
        all_blocks_clean = all_blocks;
        clean_event = event;
        return;
    end
    
    
    % start
    reject_ind = [];
    
    for i = 1:length(all_blocks.trial)
        aa = find(cell2mat(cellfun(@isnan, all_blocks.trial(i),'UniformOutput',false)));
        if ~isempty(aa)
            reject_ind(i)       = 0; % mark this trial as reject
        elseif isempty(aa)
            reject_ind(i)       = i; % don't reject this trial, keep the trial index
        end %if
    end %for

    clean_trials                = nonzeros (reject_ind)';

    % reject the marked trials
    cfg                         = [];
    cfg.trials                  = clean_trials;
    all_blocks_clean            = ft_redefinetrial(cfg, all_blocks);
    
    % also remove the indices of these rejected trials from the events
    % list, so that each index will still correspond to the correct trial!
    clean_event                 = event(clean_trials);
end