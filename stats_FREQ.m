%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_FREQ.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Grand average & statistical analysis of frequency results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% = Settings =

% Please specify correctly:
%run_name = 'offlineHPF_LMref';
run_name = 'EC_LPF30';
subj_group = 'controls'; %'migraineurs';


%%
% run the #define section
%global ResultsFolder; 
%common();

ResultsFolder = ['Z:\Analysis\Judy\EpisodicMigraine\results\' subj_group '\'];

ResultsFolder_thisrun = [ResultsFolder run_name '\\']; % where to read in the results for all subjects
save_location = [ResultsFolder_thisrun 'Figures_GA\\']; % where to save the GA figures
mkdir(save_location);

allSubjects_freq = {};


%% Read in data & compute GA (i.e. average across subjects)

%{
temp = load([ResultsFolder_thisrun '552_S1_offline0.01HPF_noICA_carefulReject.mat']);
freq1 = temp.freq;
temp = load([ResultsFolder_thisrun '9009-test_offlineHPF_noICA_carefulReject.mat']);
freq2 = temp.freq;
temp = load([ResultsFolder_thisrun '9002_S1_EC.mat']);
freq3 = temp.freq;
%}

%files = dir([ResultsFolder_thisrun '*_S1.mat']);
files = dir([ResultsFolder_thisrun 'Subject_*.mat']);

% each cycle reads in one '.mat' file (i.e. one subject's freq results)
for i = 1:length(files)
    filename = [ResultsFolder_thisrun files(i).name];
    load(filename);
    
    allSubjects_freq = [allSubjects_freq freq];
end

% export the indi-subject freq results to Excel sheet
Excel_output_file = [save_location 'summary.xls'];
if (exist(Excel_output_file, 'file') ~= 2)     
    for i = 1:length(files)
        % write the heading (SubjectID + channel labels)
        filename = [ResultsFolder_thisrun files(i).name];
        SubjectID = ['ID ' filename(end-6:end-4)];
        writecell([SubjectID 'Freq' allSubjects_freq{i}.label'], Excel_output_file, 'WriteMode','append');
        
        % write the power spectrum matrix for this subject
        M = allSubjects_freq{i}.powspctrm';
        freq_labels = 1:30;
        M2 = [NaN(height(M),1) freq_labels' M]; % add an empty col in front, then a second col containing the freq labels (1-30Hz)
        writematrix(M2, Excel_output_file, 'WriteMode','append');
    end
end

% find the channels that all subjects have
%common_chans = intersect(freq1.label, freq2.label);
%common_chans = intersect(common_chans, freq3.label);

% compute GA
cfg = [];
%cfg.foilim = [0 30];
%cfg.channel = common_chans;
%cfg.nanmean = 'yes'; % this is no use - it will only GA the common channels,
% coz we didn't keep the rejected chans as 'nan' (as that would cause problems in ft_freqanalysis)
GA_freq = ft_freqgrandaverage(cfg, allSubjects_freq{:});

% save the GA
GA_output_file = [save_location 'GA_avg.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_freq');
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


% freq3 doesn't contain infraslow (used online filter 0.3Hz), so redo here
%{
[GA_freq] = ft_freqgrandaverage([], freq1, freq2);
plot_TFR_topo(GA_freq, lay, 'infraslow', [0.03 0.06], save_location)
%}


%%

%TODO% next step is stats - see Flavia paper, no need to do clusters
        


% If you want to find spatial cluster (following FT tutorial), 
% need to check & make sure the channel layout is correct, 
% otherwise the neighbours will be incorrect


