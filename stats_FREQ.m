%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_FREQ.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Collate all subjects' results into a spreadsheet.
% Compute grand average & produce plots.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global SubjectIDs_GA; global ResultsFolder_thisrun; global LAYOUT_FILE; 
global ANALYSE_ISO; global PLOT_XLIM; global FREQ_FIELD; global FREQS_TO_EXPORT; 
global is_conn;
common();

% location to save the GA & figures
save_location = [ResultsFolder_thisrun 'GA\'];
mkdir(save_location);


%% Read in each subject's result file & plot overall power (i.e. avg of all sensors), 

% if we are working with conn results, skip this plot
if ~is_conn
    
    % plot "overall power" (one line == one subject)
    figure; hold on;
    
    % each cycle reads in one '.mat' file (i.e. one subject's freq results)
    for i = 1:length(SubjectIDs_GA)
        %filename = [ResultsFolder_thisrun files(i).name];
        filename = [ResultsFolder_thisrun cell2mat(SubjectIDs_GA(i)) '.mat'];

        load(filename);

        freq.freq = FREQ_FIELD;  % do some fixing up (coz the "freq" field created by ft_freqanalysis 
        % does not contain whole numbers & vary across subjects, causing issues for the plot)

        % plot the "overall power" for this subject
        plot(freq.freq, mean(freq.powspctrm));
        xlim(PLOT_XLIM);
        xlabel('Frequency (Hz)');
        ylabel('Absolute power (uV^2)');

        % plot the "overall power" for this subject (log transformed)
        %plot(freq.freq, mean(log(freq.powspctrm)));
        %xlim(PLOT_XLIM);
        %xlabel('Frequency (Hz)');
        %ylabel('Power (log[uV^2]');
    end
    hold off;

    % save the plot
    export_fig(gcf, [save_location 'overall_power_' int2str(length(SubjectIDs_GA)) 'subj.png']); % use this tool to save the figure exactly as shown on screen
end


%% Collate results for all subjects into one cell array

allSubjects_freq = {};

% each cycle reads in one '.mat' file (i.e. one subject's freq results)
for i = 1:length(SubjectIDs_GA)
    filename = [ResultsFolder_thisrun cell2mat(SubjectIDs_GA(i)) '.mat'];
    load(filename);
    
    % add to cell array
    if is_conn
        coh.dimord = 'subj_chan_freq'; % make a fake dimord, coz ft_freqgrandaverage
                                       % does not accept 'chan_chan_freq'
        allSubjects_freq = [allSubjects_freq coh];
    else
        allSubjects_freq = [allSubjects_freq freq];
    end
end

% save the var
save([save_location 'allSubjects_freq.mat'], 'allSubjects_freq');


%% Export the indi-subject freq results to Excel sheet

% if we are working with conn results, skip this step
% coz the 3D coherence matrix (chan x chan x freq) is difficult to export properly
if ~is_conn
    Excel_output_file = [save_location 'summary.xls'];
    if (exist(Excel_output_file, 'file') ~= 2)     
        for i = 1:length(SubjectIDs_GA)
            % write the heading (SubjectID + channel labels)
            %filename = [ResultsFolder_thisrun cell2mat(SubjectIDs(i)) '.mat'];
            %SubjectID = ['ID ' filename(end-6:end-4)];
            SubjectID = ['ID ' cell2mat(SubjectIDs_GA(i))];
            writecell([SubjectID 'Freq' allSubjects_freq{i}.label'], Excel_output_file, 'WriteMode','append');

            % write the power spectrum matrix for this subject
            M = allSubjects_freq{i}.powspctrm';
            rows = find(ismembertol(FREQ_FIELD, FREQS_TO_EXPORT, 1e-4)); % find the rows to export (comparing floating point numbers is tricky, so we add a small tolerance to make sure they can be matched up)
            M = M(rows, :);
            M2 = [NaN(height(M),1) FREQS_TO_EXPORT' M]; % add an empty col in front, then a second col containing the freq labels (e.g. 1-30Hz)
            writematrix(M2, Excel_output_file, 'WriteMode','append');
        end
    end
end


%% Compute GA (i.e. average across subjects)

% find the channels that all subjects have
%common_chans = intersect(freq1.label, freq2.label);
%common_chans = intersect(common_chans, freq3.label);

cfg = [];
%cfg.foilim = [0 30];
%cfg.channel = common_chans;
%cfg.nanmean = 'yes'; % this is no use - it will only GA the common channels,
% coz we didn't keep the rejected chans as 'nan' (as that would cause problems in ft_freqanalysis)
if is_conn
    cfg.parameter = 'cohspctrm';
end
cfg.keepindividual = 'no'; % average across subjects
GA_freq = ft_freqgrandaverage(cfg, allSubjects_freq{:});
cfg.keepindividual = 'yes'; % do not average across subjects, keep the data for each individual subject
GA_freq_indi = ft_freqgrandaverage(cfg, allSubjects_freq{:});

if is_conn
    % for conn results: we made fake dimord above, restore it to the correct info now
    GA_freq.dimord = 'chan_chan_freq';
    GA_freq_indi.dimord = 'subj_chan_chan_freq';
end

% save the GA files
GA_output_file = [save_location 'GA_avg.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_freq');
end
GA_output_file = [save_location 'GA_individuals.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_freq_indi');
end


%% GA Plots

if is_conn % for connectivity analysis
    
    % GA coherence
    figure;
    cfg           = [];
    cfg.parameter = 'cohspctrm';
    cfg.xlim      = PLOT_XLIM;
    cfg.zlim      = [0 1];
    ft_connectivityplot(cfg, GA_freq);

    set(gcf, 'Position', get(0, 'Screensize')); % make the figure full-screen
    export_fig(gcf, [save_location 'GA_coherence.png']); % use this tool to save the figure exactly as shown on screen

    % GA connectivity for each freq band
    figure;
    a = mean(GA_freq.cohspctrm(:,:, [9:12]), 3); % avg over alpha band
    imagesc(a); colorbar; ylabel('EEG channel'); xlabel('EEG channel');
    export_fig(gcf, [save_location 'alpha.png']);

    a = mean(GA_freq.cohspctrm(:,:, [13:25]), 3); % avg over beta band
    imagesc(a); colorbar; ylabel('EEG channel'); xlabel('EEG channel');
    export_fig(gcf, [save_location 'beta.png']);

    a = mean(GA_freq.cohspctrm(:,:, [4:8]), 3); % avg over theta band
    imagesc(a); colorbar; ylabel('EEG channel'); xlabel('EEG channel');
    export_fig(gcf, [save_location 'theta.png']);
    
else % for standard freq analysis (GA power spectrum & topo for each freq band)
    
    load(LAYOUT_FILE);
    plot_TFR(GA_freq, lay, save_location, PLOT_XLIM, ANALYSE_ISO); % include topoplot for infra-slow?
    
    % For sanity check: detailed topoplots (at regular freq interval)
    %{
    save_location_detailed = [save_location 'detailed_topoplots\\'];
    mkdir(save_location_detailed);

    low_freqs = [0 0.1];
    for i = 1:10
        plot_TFR_topo(GA_freq, lay, [mat2str(low_freqs(2)) 'Hz'], low_freqs, save_location_detailed);
        low_freqs = low_freqs + 0.1;
    end
    high_freqs = [1 2];
    for i = 1:29
        plot_TFR_topo(GA_freq, lay, [mat2str(high_freqs(2)) 'Hz'], high_freqs, save_location_detailed);
        high_freqs = high_freqs + 1;
    end
    %}

end
