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
run_name = 'offlineHPF';



%%
% run the #define section
global ResultsFolder; 
common();

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

files = dir([ResultsFolder_thisrun '*_S1.mat']);

% each cycle reads in one '.mat' file (i.e. one subject's freq data)
for i = 1:length(files)
    filename = [ResultsFolder_thisrun files(i).name];
    load(filename);
    
    allSubjects_freq = [allSubjects_freq freq];
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
load('lay_AntNeuro64.mat');
plot_TFR(GA_freq, lay, save_location);

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


