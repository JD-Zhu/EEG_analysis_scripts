% === exclude trials involving beh errors (based on errorsheet from manual checking) ===
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function events_allBlocks = exclude_beh_errors(SubjectID, events_allBlocks)
    % for each subject, we have 5 separate files to read:
    % {'mixedLang_natuni'}, {'mixedLang_artuni'}, {'mixedLang_bi'}, {'singleLang_Chn'}, {'singleLang_Eng'}
    
    
    % First, go thru the 3 mixed-lang contexts
    contexts = [{'mixedLang_natuni'}, {'mixedLang_artuni'}, {'mixedLang_bi'}];
    
    for i = 1:length(contexts)  
        context_name = contexts{i};
        fprintf([context_name ':\n']);
        
        % get trial info & error info for all crit trials
        allCritTrials_table = read_errorsheet_Exp2(SubjectID, context_name);
        
        % separate the full table into 4 subtables (1 for each condition)
        chstay_rows = strcmp(allCritTrials_table.lang,'Chn') & strcmp(allCritTrials_table.t_type,'Stay');
        chswitch_rows = strcmp(allCritTrials_table.lang,'Chn') & strcmp(allCritTrials_table.t_type,'Switch');
        enstay_rows = strcmp(allCritTrials_table.lang,'Eng') & strcmp(allCritTrials_table.t_type,'Stay');
        enswitch_rows = strcmp(allCritTrials_table.lang,'Eng') & strcmp(allCritTrials_table.t_type,'Switch');
        
        clear subtables;    
        subtables.chstay = allCritTrials_table(chstay_rows,:);
        subtables.chswitch = allCritTrials_table(chswitch_rows,:);
        subtables.enstay = allCritTrials_table(enstay_rows,:);
        subtables.enswitch = allCritTrials_table(enswitch_rows,:);

        % Special provisions for A07-WG-3509 (forgot to record at the
        % beginning of the first ArtUni block, so the MEG data is missing 
        % a few trials. Need to remove these trials from the error tables)
        if (strcmp(SubjectID, 'A07-WG-3509') && strcmp(context_name, 'mixedLang_artuni'))
            % sort each subtable (just in case it's not already sorted)
            subtables.chstay = sortrows(subtables.chstay, [1 2]); 
            subtables.chswitch = sortrows(subtables.chstay, [1 2]); 
            subtables.enstay = sortrows(subtables.chstay, [1 2]); 
            subtables.enswitch = sortrows(subtables.chstay, [1 2]); 
            
            % remove the trials that are missing from the MEG recording
            subtables.chstay(1:2,:) = [];
            subtables.chswitch(1,:) = [];
            subtables.enstay(1,:) = [];
            subtables.enswitch(1,:) = [];
        end

        
        % Find the corresponding eventname for each subtable
        % (1) set the appropriate prefix according to the context name
        if strcmp(context_name, 'mixedLang_natuni')
            cond_prefix = 'Nat';
        elseif strcmp(context_name, 'mixedLang_artuni')
            cond_prefix = 'Art';
        elseif strcmp(context_name, 'mixedLang_bi')
            cond_prefix = 'Bi';
        else
            fprintf('Error in exclude_beh_errors(): context_name is invalid.\n');
        end
        
        % (2) put the table names in a list, and the cond names in a list
        % These two lists must be in the same order!
        table_name_list = [{'chstay'}, {'chswitch'}, {'enstay'}, {'enswitch'}];
        cond_name_list = [{'StayC'}, {'SwitchC'}, {'StayE'}, {'SwitchE'}];      

        % Remove error trials in each cond
        for j = 1:length(cond_name_list)
            cond_name = [cond_prefix cond_name_list{j}]; % e.g. 'NatStayC'
            
            % sort the trial list (just in case it's not already sorted)
            % sorting based on round_n then trial_n
            table_sorted = sortrows(subtables.(table_name_list{j}), [1 2]); 
            
            % check to ensure #trials are the same btwn errorsheet & MEG data
            assert(height(table_sorted) == length(events_allBlocks.(cond_name)));
            
            % find the indices of trials needing to be excluded
            remove = find(table_sorted.error); % any non-zero entry indicates exclusion
            events_allBlocks.(cond_name)(remove) = []; % remove these trials
        end
    end
    
    
    % Next, we look at the single-lang blocks
    contexts = [{'singleLang_Chn'}, {'singleLang_Eng'}];
    
    for i = 1:length(contexts)  
        context_name = contexts{i};
        fprintf([context_name ':\n']);
        
        % get trial info & error info for all crit trials
        allCritTrials_table = read_errorsheet_Exp2(SubjectID, context_name);
        
        % separate the full table into 3 subtables (1 for each condition)
        natuni_rows = strcmp(allCritTrials_table.set,'NAT_UNI_C') | strcmp(allCritTrials_table.set,'NAT_UNI_E');
        artuni_rows = strcmp(allCritTrials_table.set,'ART_UNI_1') | strcmp(allCritTrials_table.set,'ART_UNI_2');
        bi_rows = strcmp(allCritTrials_table.set,'BIVALENT');

        clear subtables;    
        subtables.natuni = allCritTrials_table(natuni_rows,:);
        subtables.artuni = allCritTrials_table(artuni_rows,:);
        subtables.bi = allCritTrials_table(bi_rows,:);

        
        % Find the corresponding eventname for each subtable
        % (1) set the appropriate suffix according to the lang
        if strcmp(allCritTrials_table.lang(1),'Chn')
            cond_suffix = 'C';
        elseif strcmp(allCritTrials_table.lang(1),'Eng')
            cond_suffix = 'E';
        else
            fprintf('Error in exclude_beh_errors(): lang is invalid.\n');
        end
        
        % (2) put the table names in a list, and the cond names in a list
        % These two lists must be in the same order!
        table_name_list = [{'natuni'}, {'artuni'}, {'bi'}];
        cond_name_list = [{'Nat'}, {'Art'}, {'Bi'}];      

        % Remove error trials in each cond
        for j = 1:length(cond_name_list)
            cond_name = [cond_name_list{j} 'Single' cond_suffix]; % e.g. 'NatSingleC'
            
            % sort the trial list (just in case it's not already sorted)
            % sorting based on round_n then trial_n
            table_sorted = sortrows(subtables.(table_name_list{j}), [1 2]); 
            
            % check to ensure #trials are the same btwn errorsheet & MEG data
            assert(height(table_sorted) == length(events_allBlocks.(cond_name)));
            
            % find the indices of trials needing to be excluded
            remove = find(table_sorted.error); % any non-zero entry indicates exclusion
            events_allBlocks.(cond_name)(remove) = []; % remove these trials
        end
    end
    
end