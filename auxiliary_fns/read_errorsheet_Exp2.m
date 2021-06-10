% for a given SubjectID:
% combine its RT output (from neurobs) & errorsheet (from manual checking)
%
% Returns:
% 1 table listing all crit trials (containing trial info & error info).
% error code 4 - beh error
% error code 9 - exclude from RT analysis (due to error on prev trial)
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
%{
Version log:
v1 - only reads in original errorsheet (crit trials only)
v2 - also reads in additional errorsheet (containing manually checked
errors on filler trials) if it exists
Exp2 - for MEG exp2
%}

function RTs_crit_sorted = read_errorsheet_Exp2(SubjectID, context_name)
    
    RT_folder = [pwd '\\..\\..\\6_MEG-data\\Oak-Data\\' SubjectID '\\']; 
    RT_file = [RT_folder 'RTs_' context_name '.txt']; % output from neurobs
    
    errorsheet_folder = [pwd '\\..\\..\\6_MEG-data\\errorsheets\\' SubjectID '\\']; 
    errorsheet_file = [errorsheet_folder context_name '.csv']; % errorsheet from manual checking (crit trials only)
    errorsheet_file_fillers = [errorsheet_file(1:end-4) '_FILLERS.csv']; % this errorsheet may or may not exist
    
        
    % read in RT table
    RTs = readtable(RT_file, 'Delimiter','\t');
    %crit_rows = strcmp(RTs.t_type, 'NA')==0; % find all crit trials (ie. trial type is not 'NA')
    crit_rows = strcmp(RTs.t_type, 'Filler')==0; % find all crit trials (ie. trial type is not 'Filler')
    RTs_crit = RTs(crit_rows,:);
    RTs_crit_sorted = sortrows(RTs_crit, [1 2]); % sort the trials based on round_n then trial_n
    fprintf('  neurobs RT file found\n');

    % read in error table (contains crit trials only)
    errors = readtable(errorsheet_file);
    errors_sorted = sortrows(errors, [1 2]); % sort the trials based on round_n then trial_n
    fprintf('  standard errorsheet found\n');

    
    % ensure the RT table & error table contain same list of trials
    assert(height(RTs_crit_sorted) == height(errors_sorted));
    assert(isequal(RTs_crit_sorted.round_n, errors_sorted.round_n));
    assert(isequal(RTs_crit_sorted.trial_n, errors_sorted.trial_n));
    
    % combine errorsheet into the RT table
    RTs_crit_sorted.badTrial = errors_sorted.badTrial;
    RTs_crit_sorted.error = errors_sorted.error;


    % Additional 2-step processing to exclude all crit trials following an error on prev trial

    % Step 1 (automatic) - for all error trials, exclude the trial following it, 
    % except when the error trial itself is a switch trial (which means the
    % following trial lands on a filler)
    errors_rows = (RTs_crit_sorted.error ~= 0) & (strcmp(RTs_crit_sorted.t_type,'Switch') == 0);
    % loop thru all the identified errors trials
    for i = find(errors_rows)'
        if i < height(RTs_crit_sorted) % make sure this is not the last trial in the whole exp (in which case "i+1" exceeds matrix dimension)
            if (RTs_crit_sorted{i, 'round_n'} == RTs_crit_sorted{i+1, 'round_n'}) % make sure this is not the last trial in any block
                if (RTs_crit_sorted{i+1, 'error'} == 0) % if not already marked as error
                    RTs_crit_sorted{i+1, 'error'} = 9; % mark for exclusion (error code 9: exclude due to error on prev trial)
                end
            end
        end
    end       

    % Step 2 (based on manual checking of filler errors) - for all errors
    % that occur on a filler trial, exclude the trial after
    if (exist(errorsheet_file_fillers, 'file') == 2) % check whether the filler errorsheet exists
        fprintf('  filler errorsheet found\n');

        errors_fillers = readtable(errorsheet_file_fillers); % if so, read it in
        errors_fillers_sorted = sortrows(errors_fillers, [1 2]); % sort the trials based on round_n then trial_n
        errors_filler_rows = (errors_fillers_sorted.error ~= 0);
        % loop thru all errors occuring on a filler trial
        for i = find(errors_filler_rows)' 
            round_n = errors_fillers_sorted{i, 'round_n'}; % firstly identify this filler trial
            trial_n = errors_fillers_sorted{i, 'trial_n'}; % (i.e. find its round_n & trial_n)
            trial_n = trial_n + 1; % then identify the crit trial after (this is the trial we want to mark as exclude)
            % look for this crit trial in the big table, see which row it is on
            row = find((RTs_crit_sorted.round_n == round_n) & (RTs_crit_sorted.trial_n == trial_n));
            if (RTs_crit_sorted{row, 'error'} == 0) % if not already marked as error
                RTs_crit_sorted{row, 'error'} = 9; % mark for exclusion (error code 9: exclude due to error on prev trial)
            end
        end 
    end

%{    
    % special provision for SubjectID = M12-TH-2754
    % for whom I forgot to record MEG data in 1st block
    % so in the errorsheet, we also need to remove all trials in 1st block
    if strcmp(SubjectID, 'M12-TH-2754')
        remove_row = RTs_crit_sorted.round_n == 1; % all trials in 1st block
        RTs_crit_sorted(remove_row,:) = [];
    end
    
    % special provision for SubjectID = M09-BL-2731
    % for whom I forgot to record MEG data at the very beginning of B1
    % (resulted in the first crit trial missing)
    % so remove this trial in the errorsheet as well
    if strcmp(SubjectID, 'M09-BL-2731')
        remove_row = RTs_crit_sorted.round_n == 1 & RTs_crit_sorted.trial_n == 2;
        RTs_crit_sorted(remove_row,:) = [];
    end
%}    
    
    % convert error column to logicals (0 = correct, 1 = exclude)
    %RTs_crit_sorted.error = logical(RTs_crit_sorted.error); 

    %fprintf('read_errorsheet_Exp2 completed\n\n');
    
end