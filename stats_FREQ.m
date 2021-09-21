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
which_project = 'migraine'; % Options: 'SCI', 'migraine'

subj_group = 'controls'; % Options: 'migraineurs', 'controls'
run_name = 'EC_LPF30';
%run_name = 'offlineHPF_LMref';

% run the #define section
%global ResultsFolder; common();
ResultsFolder = ['Z:\Analysis\Judy\EpisodicMigraine\results\' subj_group '\'];

% can specify a subset of subjects to use,
% or leave empty (to use all subjs in the folder)
%SubjectIDs = {'Subject_583', 'Subject_673', 'Subject_680'};
SubjectIDs = {};


% automatic setup
ResultsFolder_thisrun = [ResultsFolder run_name '\']; % where to read in the result files for all subjects
save_location = [ResultsFolder_thisrun 'GA_' subj_group '\']; % where to save the GA & figures
mkdir(save_location);

% if subject list is empty, then use all results files in the folder
if isempty(SubjectIDs)
    SubjectIDs = dir([ResultsFolder_thisrun 'Subject_*.mat']);
    SubjectIDs = {SubjectIDs.name}; % extract the names into a cell array
end

allSubjects_freq = {};


% Read in each subject's result file & plot overall power, then collate results for all subjects into one Excel file

% plot "overall power" (i.e. avg of all sensors) for each subject,
% putting all subjects in same plot (one line == one subject)
figure; hold on;
if strcmp(which_project, 'migraine')
    x_limits = [2 30]; % for plotting (anything below 2Hz was affected by the HPF)
    freq_field = 1:30; % for fixing up the freq field (for some reason the freq values are not whole numbers)
elseif strcmp(which_project, 'SCI')
    x_limits = [0 30];
    freq_field = 0:0.005:30;
end

% each cycle reads in one '.mat' file (i.e. one subject's freq results)
for i = 1:length(SubjectIDs)
    %filename = [ResultsFolder_thisrun files(i).name];
    filename = [ResultsFolder_thisrun cell2mat(SubjectIDs(i))];

    load(filename);
   
    freq.freq = freq_field;  % do some fixing up
    
    % add to cell array
    allSubjects_freq = [allSubjects_freq freq];
    
    % plot the "overall power" for this subject
    plot(freq.freq, mean(freq.powspctrm));
    xlim(x_limits);
    xlabel('Frequency (Hz)');
    ylabel('Absolute power (uV^2)');
end
hold off;

% save the plot
%export_fig(gcf, [save_location 'overall_power_' int2str(length(SubjectIDs)) 'subj.png']); % use this tool to save the figure exactly as shown on screen

% save the var
save([save_location 'allSubjects_freq.mat'], 'allSubjects_freq');


%% Export the indi-subject freq results to Excel sheet
Excel_output_file = [save_location 'summary.xls'];
if (exist(Excel_output_file, 'file') ~= 2)     
    for i = 1:length(SubjectIDs)
        % write the heading (SubjectID + channel labels)
        filename = [ResultsFolder_thisrun cell2mat(SubjectIDs(i))];
        SubjectID = ['ID ' filename(end-6:end-4)];
        writecell([SubjectID 'Freq' allSubjects_freq{i}.label'], Excel_output_file, 'WriteMode','append');
        
        % write the power spectrum matrix for this subject
        M = allSubjects_freq{i}.powspctrm';
        freq_labels = 1:30;
        M2 = [NaN(height(M),1) freq_labels' M]; % add an empty col in front, then a second col containing the freq labels (1-30Hz)
        writematrix(M2, Excel_output_file, 'WriteMode','append');
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
cfg.keepindividual = 'no'; % average across subjects
GA_freq = ft_freqgrandaverage(cfg, allSubjects_freq{:});
cfg.keepindividual = 'yes'; % do not average across subjects, keep the data for each individual subject
GA_freq_indi = ft_freqgrandaverage(cfg, allSubjects_freq{:});

% save the GA files
GA_output_file = [save_location 'GA_avg.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_freq');
end
GA_output_file = [save_location 'GA_individuals.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_freq_indi');
end

% GA Plots (power spectrum & topo for each freq band)
%load('lay_AntNeuro64.mat');
%plot_TFR(GA_freq, lay, save_location, [1 30], true); % include topoplot for infra-slow
%load('lay_NeuroPrax32.mat');
plot_TFR(GA_freq, lay, save_location, [2 30], false);

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
