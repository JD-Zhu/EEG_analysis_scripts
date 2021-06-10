    % loop thru each ROI
    ROIs_label = fieldnames(ROI_activity);
    for k = 1:length(ROIs_label)
        ROI_name = ROIs_label{k};
        
        % loop thru each cond
        conds_label = fieldnames(ROI_activity.(ROI_name));
        for j = 1:length(conds_label)
            cond_name = conds_label{j};
            
            % read the original (ie. not yet smoothed) timecourse
            x = ROI_activity.(ROI_name).(cond_name).time;
            y = ROI_activity.(ROI_name).(cond_name).avg;
            
            % apply various smoothing filters & plot the resulting waveforms
            %
            figure; hold on;
            plot(x, y); % plot the original
            
            % median filter
            %
            %plot(x, medfilt1(y, 2)); % very similar to original waveform
            %plot(x, medfilt1(y, 9)); % https://www.ncbi.nlm.nih.gov/pubmed/32024267
            plot(x, medfilt1(y, 3)); % adjusted window sizes based on sampling rate (theirs is 512Hz, mine is downsampled to 200Hz)
            plot(x, medfilt1(y, 4));
            %plot(x, smoothdata(y, 'movmedian')); % can compute an appropriate "order" (ie. window size) automatically.
                                                 % like the options above, this also creates weird "horizontal stretches" in resulting waveform: 
                                                 % sometimes similar to (or better than) order 6, sometimes worse than order 9
            legend({'orig', 'order3', 'order4'});
            %}
            
            % basic smooth filter & sgolay filter
            %
            plot(x, smooth(y)); % basic smooth filter
            plot(x, sgolayfilt(y, 3, 11)); % https://www.mathworks.com/help/signal/ref/sgolayfilt.html
            %plot(x, sgolayfilt(y, 4, 27)); % https://www.ncbi.nlm.nih.gov/pubmed/32024267
            plot(x, sgolayfilt(y, 4, 11)); % adjusted window size based on sampling rate (theirs is 512Hz, mine is downsampled to 200Hz)
            legend({'orig', 'basic-smooth', 'sgolayfilt-example', 'sgolayfilt-paper'});
            %}
            hold off;
        end
    end
    