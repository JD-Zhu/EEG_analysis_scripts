% Identify the trial numbers belonging to each condition (ie. event type), 
% as defined in the eventcodes{} array

function events_allBlocks = identify_event_types(SubjectID, trialinfo_b)
    
    % special provision for SubjectID = B12-CD-3388
    % for whom I forgot to record MEG data at the beginning of B2
    % (the first response may have been partially cut off, so just remove it)
    if strcmp(SubjectID, 'B12-CD-3388')
        trialinfo_b.trl(729,:) = []; % identified this trial number manually
        trialinfo_b.event(729) = [];
    end
    

    % run the #define section
    global eventcodes; global eventnames;
    common;
    
    trialsgone = 0;

    % each cycle handles one "block" (i.e. one '.con' file)
    for i = 1:length(trialinfo_b)
        for j = 1:length(eventnames)
            events.(eventnames{j}) = find(strcmp({trialinfo_b.event.value}, eventcodes{j,2})); 
            % 9 fields representing the 9 types of events
            % each field contains a list of all events belonging to this type (by matching event code)
            % NB: this is a tmp var, it gets overwritten in each cycle
        end

        if i == 1 % first block

            for j = 1:length(eventnames)
                % save the lists to a perm var, also transpose each list
                events_allBlocks.(eventnames{j}) = events.(eventnames{j})'; 
            end

        else % all other blocks
            trialsinblock = length(trialinfo_b(i-1).event); % how many "trials" (i.e. events) were identified in previous block
            trialsgone = trialsgone + trialsinblock; % add this number to the total number of "past" trials

            for j = 1:length(eventnames)
                events_allBlocks.(eventnames{j}) = [events_allBlocks.(eventnames{j}); events.(eventnames{j})' + trialsgone]; 
                % continue to append to the perm lists (stored in "events_allBlocks")                                                                                             
                % in the end, each perm list will contain all events of that type from all blocks
            end
        end
    end

    % special provision for SubjectID = M09-BL-2731
    % for whom I forgot to record MEG data at the very beginning of B1
    % (the first response may have been partially cut off, so just remove it)
    if strcmp(SubjectID, 'M09-BL-2731')
        events_allBlocks.response(1) = [];
    end

end