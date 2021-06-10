time = {trialinfo.event.time}'; % a list of time of all events (cue, target, photodetector, response)
time_event = time(find(ismember({trialinfo.event.value}','5'))-1); % all target events (ie. the event b4 each photodetector)
time_photo = time(find(ismember({trialinfo.event.value}','5'))); % all photodetector events
time_diff = cell2mat(time_photo)-cell2mat(time_event);
figure; plot(time_diff,'.'); %plot the diff btwn target & (the immediately following turning-on of) photodetector