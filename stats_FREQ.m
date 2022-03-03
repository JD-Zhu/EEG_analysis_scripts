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


% = Settings =

% PLEASE SPECIFY:
which_project = 'SCI'; % Options: 'SCI', 'migraine'
subj_group = ''; % Options: 'migraineurs', 'controls' (for SCI project, please leave empty for now)

%ProjectFolder = 'Z:\Analysis\Judy\EpisodicMigraine\';
ProjectFolder = 'Z:\Analysis\Preprocess\NeuRA_SCI_SCS_CIPN_BUMP\EEG\';
%run_name = 'EC_LPF30';
run_name = 'offlineHPF_LMref';

% are we working with connectivity results here?
is_conn = false;


% can specify a subset of subjects to use,
% or leave empty (to use all subjs in the folder)
SubjectIDs = [];
% Final set of 17 controls (age & gender matched to migraineurs)
if strcmp(which_project, 'migraine') && strcmp(subj_group, 'controls')
    SubjectIDs = {'Subject_101', 'Subject_251', 'Subject_252', 'Subject_253', 'Subject_254', 'Subject_495', 'Subject_610', 'Subject_622', 'Subject_623', 'Subject_634', 'Subject_642', 'Subject_675', 'Subject_690', 'Subject_809', 'Subject_844', 'Subject_885', 'Subject_891'};
end
% Groups based on migraine phases:
%SubjectIDs = {'Subject_500', 'Subject_548', 'Subject_208'}; % prodrome
%SubjectIDs = {'Subject_583', 'Subject_673', 'Subject_680', 'Subject_205'}; % postdrome
%SubjectIDs = {'Subject_661', 'Subject_664', 'Subject_671', 'Subject_677', 'Subject_681', 'Subject_696', 'Subject_800', 'Subject_207', 'Subject_209', 'Subject_210'}; % interictal
% Groups based on migraine frequency:
%SubjectIDs = {'Subject_677', 'Subject_681', 'Subject_696', 'Subject_800'}; % <1 day / month
%SubjectIDs = {'Subject_583', 'Subject_661', 'Subject_671'}; % 1-2 days / month
%SubjectIDs = {'Subject_500', 'Subject_548', 'Subject_664', 'Subject_673', 'Subject_680'}; % >3 days / month


% automatic setup
%global ResultsFolder; common();
ResultsFolder = [ProjectFolder 'results\' subj_group '\']; % all subjects' freq analysis results are stored here
if is_conn
    ResultsFolder = [ProjectFolder 'results_conn\' subj_group '\']; % all subjects' connectivity results are stored here
    %run_name = [run_name '_afterSL']; % if you want to use the version of results with SL applied
end
ResultsFolder_thisrun = [ResultsFolder run_name '\']; % where to read in the result files for all subjects
save_location = [ResultsFolder_thisrun 'GA_' subj_group '\']; % where to save the GA & figures
mkdir(save_location);

% if subject list is empty, then use all results files in the folder
if isempty(SubjectIDs)
    SubjectIDs = dir([ResultsFolder_thisrun '*.mat']);
    SubjectIDs = {SubjectIDs.name}; % extract the names into a cell array
    SubjectIDs = cellfun(@(x) x(1:end-4), SubjectIDs, 'un', 0); % remove the '.mat' extension
end

% settings for each project
if strcmp(which_project, 'migraine')
    x_limits = [2 30]; % for plotting, we are interested in 2-30Hz (everything else was filtered out)
    freq_field = 1:30; % for fixing up the freq field (for some reason the freq values are not whole numbers)
elseif strcmp(which_project, 'SCI')
    x_limits = [1 30]; % anything below 1Hz is way over powered (rendering the whole plot unviewable)
    freq_field = 0:0.01:30;
end


%% Read in each subject's result file & plot overall power (i.e. avg of all sensors), 

% if we are working with conn results, skip this plot
if ~is_conn
    
    % plot "overall power" (one line == one subject)
    figure; hold on;
    
    % each cycle reads in one '.mat' file (i.e. one subject's freq results)
    for i = 1:length(SubjectIDs)
        %filename = [ResultsFolder_thisrun files(i).name];
        filename = [ResultsFolder_thisrun cell2mat(SubjectIDs(i)) '.mat'];

        load(filename);

        freq.freq = freq_field;  % do some fixing up (coz the "freq" field created by ft_freqanalysis 
        % does not contain whole numbers & vary across subjects, causing issues for the plot)

        % plot the "overall power" for this subject
        plot(freq.freq, mean(freq.powspctrm));
        xlim(x_limits);
        xlabel('Frequency (Hz)');
        ylabel('Absolute power (uV^2)');

        % plot the "overall power" for this subject (log transformed)
        %plot(freq.freq, mean(log(freq.powspctrm)));
        %xlim(x_limits);
        %xlabel('Frequency (Hz)');
        %ylabel('Power (log[uV^2]');
    end
    hold off;

    % save the plot
    export_fig(gcf, [save_location 'overall_power_' int2str(length(SubjectIDs)) 'subj.png']); % use this tool to save the figure exactly as shown on screen
end


%% Collate results for all subjects into one cell array

allSubjects_freq = {};

% each cycle reads in one '.mat' file (i.e. one subject's freq results)
for i = 1:length(SubjectIDs)
    filename = [ResultsFolder_thisrun cell2mat(SubjectIDs(i)) '.mat'];
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
        for i = 1:length(SubjectIDs)
            % write the heading (SubjectID + channel labels)
            %filename = [ResultsFolder_thisrun cell2mat(SubjectIDs(i)) '.mat'];
            %SubjectID = ['ID ' filename(end-6:end-4)];
            SubjectID = ['ID ' cell2mat(SubjectIDs(i))];
            writecell([SubjectID 'Freq' allSubjects_freq{i}.label'], Excel_output_file, 'WriteMode','append');

            % write the power spectrum matrix for this subject
            M = allSubjects_freq{i}.powspctrm';
            freq_labels = freq_field;
            if strcmp(which_project, 'SCI') % for SCI proj, only export certain freqs (coz we computed 6001 freq points: 0:0.005:30)
                freq_labels = [0.1:0.1:0.9 1:30];
                rows = find(ismembertol(freq_field, freq_labels, 1e-4)); % find the rows to export (comparing floating point numbers is tricky, so we add a small tolerance to make sure they can be matched up)
                M = M(rows, :);
            end
            M2 = [NaN(height(M),1) freq_labels' M]; % add an empty col in front, then a second col containing the freq labels (e.g. 1-30Hz)
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
    cfg.xlim      = x_limits;
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
    
    if strcmp(which_project, 'SCI')
        load('lay_AntNeuro64.mat');
        plot_TFR(GA_freq, lay, save_location, x_limits, true); % include topoplot for infra-slow
    elseif strcmp(which_project, 'migraine')
        load('lay_NeuroPrax32.mat');
        plot_TFR(GA_freq, lay, save_location, x_limits, false);
    end
    
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
