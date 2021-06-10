%{
Version log:

v2: photodetector channel reads triggers both going up & down;
    first 5 subjects code_delay = 30ms, all other subjects code_delay = 100ms

%}

function [trl, event] = trig_fun_160_basic_v2(cfg)
    global SubjectID; % access SubjectID from base workspace
    global ADJUST_SCREEN_DELAY;
    common;
    
    
    % Change as required:
    PHOTODETECTOR_CHANNEL = 48; % 63;  % old: 5
    TRIGGER_CHANNELS = 17:95; 
    VOICEKEY_CHANNEL = 30;
    
    %{
    subject = SubjectID(1:6);
    if (strcmp(subject, 'M04-LL') || strcmp(subject, 'M05-XY') || ...
        strcmp(subject, 'M07-YH') || strcmp(subject, 'M09-BL') || strcmp(subject, 'M11-LL')) % these 5 subjects have code_delay set to 30ms
        code_delay = 0.03; % 30ms
    else % all newer subjects have code_delay set to 100ms (due to the "ignore short response" setting)
        code_delay = 0.1; % 100ms
    end
    fprintf('SubjectID = %s\ncode_delay = %d ms\n', SubjectID, code_delay*1000);
    %}
    CODE_DELAY = 0.1; % 100ms (we chose this value in Presentation)
    
    
    %% start
    all_trigger_channels  =  161:225; % channels 160-224 are used for triggers
    triggers          = ft_read_data(cfg.dataset, 'dataformat','yokogawa_con', 'chanindx',all_trigger_channels); % read from the trigger channels
    triggers(33,:)    = 0; % channel 192 is normally high, and goes down at all events (not used for trigger here)
    %triggers(7,:)    = 0; % channel 166 is the audio channel (not used for trigger here)
    triggers(64,:)    = 0; % channel 223 is the audio channel (not used for trigger here)
    triggers(1:4,:)   = 0; % these other channels are also not used for trigger here,
    triggers(8:16,:)  = 0; % some of them (e.g. the ref sensors) are crossing thresh, 
                         % so we just remove these channels to avoid them being treated as trigger events
    if strcmp(SubjectID, 'B01-YS-3479') % only for Yanan's pilot data (lots of flickering on photodetector channel, 
        triggers(64,:)   = 0;           % so just remove altogether & use the default delay)
    end
        
        
    hdr = ft_read_header(cfg.dataset, 'dataformat','yokogawa_con');
    
    for i=1:size(triggers,1)
        trig_height(i) = max(triggers(i,:)); %Y = prctile(x,42)
    end
    
    trig_thresh = 0.25 * max(trig_height); % set threshold to 1/4 of max trigger height (5V)
    fprintf('trigger threshold is %d Volt\n', trig_thresh);
    triggers    = triggers > trig_thresh; %Binarise trigger channels (1 for above threshold; 0 for no event)
     
    event = []; % this will contain a single list of all events (from all channels)
    for i = 1:size(triggers,1) % each cycle is one channel
        channel   = all_trigger_channels(i);
        trig      = triggers(i,:); % all the 1s and 0s for this channel
        pad       = trig(1);
        trigshift = 0; % SHOULD THIS BE SET TO THE PULSE_WIDTH?
        begsample = 1;
                
        changes = diff([pad trig(:)']); % find all changes in this channel (i.e. flip from 0 to 1, flip from 1 to 0)
%{
       if (i == 6) % for the photodetector channel, take both +ve and -ve changes (trigger going up/down)
                    % +ve change (trigger going up) is white box, -ve (going down) is black box
            changes_index = find(changes < 0); %  ~=
        else % for all other channels, only take all the +ve changes (trigger going up)
            changes_index = find(changes > 0);
        end
%}
        % for all trigger channels (incl. photodetector channel), only take the +ve changes (trigger going up)
        changes_index = find(changes > 0);
         
        for j = changes_index % for each relevant trigger channel voltage change (up/down), add this as a trigger event
            event(end+1).type   = num2str(channel);       % add this "event" (channel#, sample#, value) to the end of the "event" array
            event(end  ).sample = j + begsample - 1;      % the sample at which the trigger has gone down
            event(end  ).value  = trig(j + trigshift);    % the trigger value just _after_ going up
        end    
    end
    
    %% Fill up the event struct
    event2 = nestedSortStruct(event,'sample'); % sort all trigger events by chronological order
    event3 = {event2.type}; % channel number for each event
    event3 = cell2mat(cellfun(@str2num,event3,'un',0)) - 161; % "channel numbers" now corrected to correspond to port_codes
    event3 = arrayfun(@num2str, event3, 'unif', 0);
    event4 = repmat(struct('type','trigger','value',1,'sample',1,'time',1,'duration',''), 1, length(event3));
    
    times2 = [event2.sample];
    
    for m=1:length(event3)
        event4(m).time   = times2(m)/hdr.Fs;
        event4(m).sample = times2(m);
        event4(m).value  = cell2mat(event3(m));
    end
    
    event = event4;
    
    
    %%%% TRIGGER ADJUSTMENTS FOR CODE_DELAY, PHOTODETECTOR, AND
    %%%% REASSIGN CODES BASED ON TASK
    
    events1 = event;
    events2 = {events1.value};
    
    events2     = cell2mat(cellfun(@str2num,events2,'un',0));
    event_times = [events1.time]';
    
    trigger1_idx = min(find(ismember(events2, TRIGGER_CHANNELS))); % trial type triggers (port codes 17~24)
    events2      = events2(trigger1_idx:end)';
    event_times  = event_times(trigger1_idx:end);
    
    if isempty(find(ismember(events2, PHOTODETECTOR_CHANNEL))) % no photodetector (channel 165)
        photodetectorcodes_n = 0;
        fprintf('No photodetector found\n');
    else % if there are photodetector triggers, count how many there are
        photodetectorcodes_n = sum(ismember(events2, PHOTODETECTOR_CHANNEL));
    end
fprintf('photodetectorcodes_n = %d\n',photodetectorcodes_n);
    
    if photodetectorcodes_n > 20 %isalmost(photodetectorcodes_n,gocodes_n,10)  % if there are at least 20 photodetector triggers, then use it to adjust timing
        photodetector_times = event_times(ismember(events2, PHOTODETECTOR_CHANNEL)); % all photodetector trigger times
%fprintf('how many photodetector trigger times found? %d\n',length(photodetector_times));
        go_times            = event_times(find(ismember(events2, PHOTODETECTOR_CHANNEL)) - 1); % the events just before each photodetector (should always be a 'target' trigger)
        %TODO: do not include photodetector channel when looking for the
        %event just before -> this change will remove all flickering problems!!
        delays              = photodetector_times - go_times; % a list of all the screen delay values for individual triggers
fprintf('delays array has length of %d\n', length(delays));

        delays(delays > 0.2) = []; % screen delay should not be >200ms or <15ms. If there are some, these are prob being mapped to the wrong triggers 
        delays(delays < 0.015) = []; % (can happen when there are redundant photodetector triggers), so remove these
fprintf('delays array (after removing any delays >200ms or <15ms) has length of %d\n',length(delays));
%save(['..\\..\\screen_delays\\' SubjectID '_delays'], 'delays', '-v7.3');

        screen_delay        = median(delays); % in the list of time diffs btwn photodetector onset & target onset, 
                                              % find the mode (i.e. most-occuring value), adjust all targets by that value
    else  % otherwise, use a default delay value (56ms)
        screen_delay = 0.071;
    end

    fprintf('screen_delay (median) is %d ms\n', screen_delay*1000);

    % Chn stay: cue - channel 17, target - channel 18
    % Chn sw:	  cue - channel 19, target - channel 20
    % Eng stay: cue - channel 21, target - channel 22
    % Eng sw:	  cue - channel 23, target - channel 24
    
    
    % adjust for screen delay
    if strcmp(ADJUST_SCREEN_DELAY, 'median')
        % adjust trigger times based on the median screen delay
        event_times_adjusted = event_times + single(ismember(events2, setdiff(TRIGGER_CHANNELS, VOICEKEY_CHANNEL)))*screen_delay; 
    elseif strcmp(ADJUST_SCREEN_DELAY, 'individual_trial')
        % adjust trigger time on each indi trial based on the photodetector trigger on that trial
        event_times_adjusted = event_times;
        counter2 = 1; % this counter acts as index for the "photodetector_times" array
        for counter1 = 1:length(event_times_adjusted) % this counter acts as index for the "event_times_adjusted" array
            if ismember(events2(counter1), setdiff(TRIGGER_CHANNELS, VOICEKEY_CHANNEL)) % if this event is an actual trial (i.e. trigger channels 17:95, excluding the voice key channel 30)
                event_times_adjusted(counter1) = photodetector_times(counter2); % then set the time to be the same as the photodetector trigger time
                counter2 = counter2 + 1; % manually increment this counter
            end
        end
        % sanity check
        assert (counter2 == length(photodetector_times) + 1);
    else % should never be here
        fprintf('Error in trig_fun_160_basic_v2(): global var ADJUST_SCREEN_DELAY is set to %s, this is an invalid option.\n', ADJUST_SCREEN_DELAY);
    end
    
    %event_times_adjusted = event_times_adjusted-single(ismember(events2, 17:30))*0.75; % go back to cue onset (750ms earlier)
    
    % adjust for code delay (which affects the voice key trigger - on port code 30)
    event_times_adjusted = event_times_adjusted - single(ismember(events2, VOICEKEY_CHANNEL))*CODE_DELAY; 
    
    
    events2 = arrayfun(@num2str, events2, 'unif', 0);
    
    events3 = repmat(struct('type','trigger','value',1,'time',1,'sample',1,'duration',''), 1, length(events2));
    
    for m=1:length(events2)
        events3(m).time   = event_times_adjusted(m);
        events3(m).sample = event_times_adjusted(m)*hdr.Fs;
        events3(m).value  = cell2mat(events2(m));
        events3(m).type   = events3(m).type;
    end
    
    event = events3;
    
    % now we remove the photodetector events from the list, we only want to
    % return the "real events" (cue, target, response)
    value_column = {event.value}; % value indicates the trigger channel
    real_event_rows = str2double(value_column) ~= PHOTODETECTOR_CHANNEL; % do not include the photodectector events
    event = event(real_event_rows); % only select the rows containing "real events"  
    
    % search for "trigger" events
    value  = {event(find(strcmp('trigger', {event.type}))).value}';
    sample = round([event(find(strcmp('trigger', {event.type}))).sample]');
    
    % determine the number of samples before and after the trigger
    pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
    posttrig = round(cfg.trialdef.poststim * hdr.Fs);
    
    trl = [];
    for n = 1:length(value)
        trg1     = value(n);
        trlbegin = sample(n) + pretrig;
        trlend   = sample(n) + posttrig;
        offset   = pretrig;
        newtrl   = [trlbegin trlend offset];
        trl      = [trl; newtrl];
    end
    
end

