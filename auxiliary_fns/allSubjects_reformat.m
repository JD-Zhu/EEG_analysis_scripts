% For the given ROI, reformat the ROI activity stored in per-subject organisation 
% into a by-condition organistion.

% INPUT: a list of ROI_activity structs, one struct for each subject
% OUTPUT: a list of conditions, each containing all subjects' ROI activity
% in that condition

% This function only deals with one ROI at a time (specified by the input param).
% To process another ROI, call this function again.
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function allSubjects_ROIs_byCondition = allSubjects_reformat(allSubjects_ROIs_bySubjects, ROI_name, eventnames_real)
    % each cycle processes one condition
    for j = 1:length(eventnames_real)
        allSubjects_ROIs_byCondition.(eventnames_real{j}) = {}; % initialise the array of subjects for this cond
        % loop thru all subjects & extract their ROI activity in this condition
        for i = 1:length(allSubjects_ROIs_bySubjects)
           allSubjects_ROIs_byCondition.(eventnames_real{j}) = [allSubjects_ROIs_byCondition.(eventnames_real{j}) allSubjects_ROIs_bySubjects(i).(ROI_name).(eventnames_real{j})];
        end
    end    
end