% Computes the margin for "shaded boundary" in timecourse plots
% A margin value is calc'd at every time point
%
% Each call to this fn handles one cond (i.e. calcs margin for a single timecourse line)
%
% @param margin_type can be: STDEV, SEM, CI_95
% @param allsubjects: the GA (with keepindividual) for this condition
%
% Author: Judy Zhu (github.com/JD-Zhu)
% Based on code written by Yanan Sun
%
function margin = calc_margin(allsubjects, margin_type)
    % standard deviation
    SD = std(allsubjects);

    % standard error of the mean
    SEM = SD ./ sqrt(size(allsubjects, 1));
    %sem = squeeze(sem(1,:,:));

    % 95% CI
    CI_95 = SEM * 1.96;

    % check settings: which margin_type did we choose?
    if strcmp(margin_type, 'STDEV')
        margin = SD;
    elseif strcmp(margin_type, 'SEM')
        margin = SEM;
    elseif strcmp(margin_type, 'CI_95')
        margin = CI_95;
    else
        margin = [];
        disp('\nError: "margin_type" was incorrectly specified. Available options are: STDEV, SEM, CI_95\n');
    end
end