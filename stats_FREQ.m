%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_FREQ.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Grand average & statistical analysis of frequency results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ResultsFolder; 
common();


% Please specify correctly:
run_name = 'offlineHPF';

ResultsFolder_thisrun = [ResultsFolder run_name '\\']; % results for all subjects


%% Plot grand ave across subjects

dir([ResultsFolder_thisrun '*.mat']);

temp = load([ResultsFolder_thisrun '552_S1_offline0.01HPF_noICA_carefulReject.mat']);
freq1 = temp.freq;
temp = load([ResultsFolder_thisrun '9009-test_offlineHPF_noICA_carefulReject.mat']);
freq2 = temp.freq;
temp = load([ResultsFolder_thisrun '9002_S1_EC.mat']);
freq3 = temp.freq;

% find the channels that all subjects have
%common_chans = intersect(freq1.label, freq2.label);
%common_chans = intersect(common_chans, freq3.label);

cfg = [];
%cfg.foilim = [0 30];
%cfg.channel = common_chans;
%cfg.nanmean = 'yes'; % this is no use - it will only GA the common channels,
% coz we didn't keep the rejected chans as 'nan' (which causes problems in ft_freqanalysis)
[grandavg] = ft_freqgrandaverage(cfg, freq1, freq2, freq3);
    
% where to save the figures
save_location = [ResultsFolder_thisrun 'Figures_GA\\' run_name '\\'];
mkdir(save_location);

% plot power spectrum for all channels (overlay)
figure; hold on;
for chan = 1:length(grandavg.label)
    plot(grandavg.freq, grandavg.powspctrm(chan,:))
end
xlim([1 30]);
xlabel('Frequency (Hz)');
ylabel('Absolute power (uV^2)');
hold off;

export_fig(gcf, [save_location 'powspctrm_allchans.png']); % use this tool to save the figure exactly as shown on screen

% plot avg of all channels (log transformed)
figure; plot(grandavg.freq, mean(log(grandavg.powspctrm)));
xlim([0.01 30]);
xlabel('Frequency (Hz)');
ylabel('Power (log[uV^2]');

export_fig(gcf, [save_location 'powspctrm_avg.png']); % use this tool to save the figure exactly as shown on screen
        
% plot topography for each freq band
plot_TFR_topo(grandavg, lay, 'infraslow', [0.03 0.06], save_location)
plot_TFR_topo(grandavg, lay, 'theta', [4 8], save_location)
plot_TFR_topo(grandavg, lay, 'alpha', [9 12], save_location)
plot_TFR_topo(grandavg, lay, 'beta', [13 25], save_location)

% freq3 doesn't contain infraslow (used online filter 0.3Hz), so redo here
[grandavg] = ft_freqgrandaverage([], freq1, freq2);
plot_TFR_topo(grandavg, lay, 'infraslow', [0.03 0.06], save_location)


%%

%TODO% next step is stats - see Flavia paper, no need to do clusters
        


% If you want to find spatial cluster (following FT tutorial), 
% need to check & make sure the channel layout is correct, 
% otherwise the neighbours will be incorrect


