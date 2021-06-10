% Computes ERFs & covariance matrices for by-condition data, combined cue & combined target
%
% @param trials_clean: a struct containing all the conds (each cond is one
%     field, which contains data for all trials belonging to this cond)
% @param CALC_COV:     whether to calculate the covariance matrix
%
function [erf_clean, erf_allconds] = compute_ERF (trials_clean, CALC_COV)
    % run the #define section
    %global conds_cue; global conds_target; 
    %global conds;
    %common();
    
    % grab the list of conds in this dataset
    conds = fieldnames(trials_clean);
    
    % We kept bad trials as NaN, so need to exclude them now
    for j = 1:length(conds)
        %length(trials_clean.(conds{j}).trial)
        good_trials_idx = find(~isnan(cell2mat(cellfun(@(isgood)isgood(1),trials_clean.(conds{j}).trial,'uni',0)))); % Just need to evaluate the first element as all samples in bad trial are NaN
        cfg             = [];
        cfg.trials      = good_trials_idx;
        trials_clean.(conds{j}) = ft_redefinetrial(cfg, trials_clean.(conds{j}));
    end

    % An alternative method to remove NaN trials - doesn't seem to work!
    %{
    for j = 1:length(conds)
        cfg         = [];
        cfg.nanmean = 'yes';
        trials_clean.(conds{j}) = ft_selectdata(cfg, trials_clean.(conds{j})); % Do this because we kept bad trials as NaN
    end
    %}


    % Compute erf & cov matrix on the combined data (all conds appended together)
    cellarray = struct2cell(trials_clean); % convert struct to cell array, then you can feed it in as 'varargin'
    trials_clean_allconds = ft_appenddata([], cellarray{:});

    cfg                  = [];
    if (CALC_COV)
        cfg.covariance       = 'yes';
        cfg.covariancewindow = [0 0.5]; % do not include any period after vocal response onset
    end
    erf_allconds  = ft_timelockanalysis(cfg, trials_clean_allconds);

    % Alternative method (not sure if this is 100% correct)
    %{
    % average ERF across all conds
    erf_allconds = erf_clean.(conds{1});
    for j = 2:length(conds) % add together the ERF from all conds
        erf_allconds.avg = erf_allconds.avg + erf_clean.(conds{j}).avg;    
    end
    erf_allconds.avg = erf_allconds.avg / length(conds); % find the mean
    %}
    
    % Compute erf & cov matrix for each condition 
    for j = 1:length(conds)
        erf_clean.(conds{j}) = ft_timelockanalysis(cfg, trials_clean.(conds{j}));
    end

end