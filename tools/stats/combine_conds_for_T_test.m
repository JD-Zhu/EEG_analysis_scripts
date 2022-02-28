% combine_conds_for_T_test.m 
%
% Reduce multiple conditions down to 2 conditions, so that they can be
% submitted to ft_timelockstatistics / ept_TFCE for statistical comparison.
%
% @param data_format: 'eeglab' or 'fieldtrip'
%    'eeglab' format = 3d matrix of values (subj x chan x time)
% 'fieldtrip' format = cell array of subjects, each cell containing the timelock struct for one subject (.avg field contains chan x time matrix)
%
% @param type_of_effect: 'main_12vs34', 'main_13vs24', or 'interaction',
% corresponding to the 3 types of effect which can be tested in a 2x2 design:
% - main effect of lang (1+2 vs 3+4)
% - main effect of ttype (1+3 vs 2+4)
% - interaction (2-1 vs 4-3)
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function [timelock1, timelock2] = combine_conds_for_T_test(data_format, type_of_effect, chstay, chswitch, enstay, enswitch)
    % what is the input data format?
    if strcmp(data_format, 'fieldtrip')
        % which type of effect do you want to test? 
        if strcmp(type_of_effect, 'interaction')
            [timelock1, timelock2] = interaction_2x2(chstay, chswitch, enstay, enswitch);
        elseif strcmp(type_of_effect, 'main_12vs34')
            [timelock1, timelock2] = main_12vs34(chstay, chswitch, enstay, enswitch);
        elseif strcmp(type_of_effect, 'main_13vs24')
            [timelock1, timelock2] = main_13vs24(chstay, chswitch, enstay, enswitch);
        end
    elseif strcmp(data_format, 'eeglab') 
        % which type of effect do you want to test? 
        if strcmp(type_of_effect, 'interaction')
            [timelock1, timelock2] = interaction_2x2_eeglab(chstay, chswitch, enstay, enswitch);
        elseif strcmp(type_of_effect, 'main_12vs34')
            [timelock1, timelock2] = main_12vs34_eeglab(chstay, chswitch, enstay, enswitch);
        elseif strcmp(type_of_effect, 'main_13vs24')
            [timelock1, timelock2] = main_13vs24_eeglab(chstay, chswitch, enstay, enswitch);
        end
    else
        fprintf('Error: in combine_conds_for_T_test: data_format not specified.');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% === For data in 'fieldtrip' format ===

% Prepare the data for testing the interaction effect in a 2x2 design
% (e.g. to test lang*sw interaction: calc sw$ in each lang, then return the 2 sw$ for comparison)
% http://www.fieldtriptoolbox.org/faq/how_can_i_test_an_interaction_effect_using_cluster-based_permutation_tests
function [ch_switchCost, en_switchCost] = interaction_2x2(chstay, chswitch, enstay, enswitch)
    % initialise the 2 sw$ structs
    ch_switchCost = chstay;
    en_switchCost = enstay;

    % calc the 2 sw$ for each subject
    numSubjects = length(chstay);
    for i = 1:numSubjects
        ch_switchCost{i}.avg = chswitch{i}.avg - chstay{i}.avg; % chswitch - chstay
        en_switchCost{i}.avg = enswitch{i}.avg - enstay{i}.avg; % enswitch - enstay
    end
end

% Prepare the data for testing the main effect of 1+2 vs 3+4
% (i.e. calc the average of 1 & 2, calc the average of 3 & 4, then compare)
function [ch, en] = main_12vs34(chstay, chswitch, enstay, enswitch)
    % initialise ch & en
    ch = chstay;
    en = enstay;
    
    % calc the average for ch & average for en
    numSubjects = length(chstay);
    for i = 1:numSubjects
        ch{i}.avg = (chstay{i}.avg + chswitch{i}.avg) / 2; % average across stay & switch
        en{i}.avg = (enstay{i}.avg + enswitch{i}.avg) / 2; % average across stay & switch
    end
end

% Prepare the data for testing the main effect of 1+3 vs 2+4
% (i.e. calc the average of 1 & 3, calc the average of 2 & 4, then compare)
function [stay, sw] = main_13vs24(chstay, chswitch, enstay, enswitch)
    % initialise stay & sw
    stay = chstay;
    sw = chswitch;
    
    % calc the average for stay & average for switch
    numSubjects = length(chstay);
    for i = 1:numSubjects
        stay{i}.avg = (chstay{i}.avg + enstay{i}.avg) / 2; % average across ch & en
        sw{i}.avg = (chswitch{i}.avg + enswitch{i}.avg) / 2; % average across ch & en
    end
end


% === For data in 'eeglab' format ===

% Prepare the data for testing the interaction effect in a 2x2 design
% (e.g. to test lang*sw interaction: calc sw$ in each lang, then return the 2 sw$ for comparison)
% http://www.fieldtriptoolbox.org/faq/how_can_i_test_an_interaction_effect_using_cluster-based_permutation_tests
function [ch_switchCost, en_switchCost] = interaction_2x2_eeglab(chstay, chswitch, enstay, enswitch)
    % calc the 2 sw$ for each subject
    ch_switchCost = chswitch - chstay;
    en_switchCost = enswitch - enstay;
    
    % same effect as the manual looping below
    %{
    numSubjects = size(chstay, 1);
    for i = 1:numSubjects
        ch_switchCost(i,:,:) = chswitch(i,:,:) - chstay(i,:,:); % chswitch - chstay
        en_switchCost(i,:,:) = enswitch(i,:,:) - enstay(i,:,:); % enswitch - enstay
    end
    %}
end

% Prepare the data for testing the main effect of 1+2 vs 3+4
% (i.e. calc the average of 1 & 2, calc the average of 3 & 4, then compare)
function [ch, en] = main_12vs34_eeglab(chstay, chswitch, enstay, enswitch)
    % calc the average for ch & average for en
    ch = (chstay + chswitch) / 2; % average across stay & switch
    en = (enstay + enswitch) / 2; % average across stay & switch
end

% Prepare the data for testing the main effect of 1+3 vs 2+4
% (i.e. calc the average of 1 & 3, calc the average of 2 & 4, then compare)
function [stay, sw] = main_13vs24_eeglab(chstay, chswitch, enstay, enswitch)
    % calc the average for stay & average for switch
    stay = (chstay + enstay) / 2; % average across ch & en
    sw = (chswitch + enswitch) / 2; % average across ch & en
end
