%% Run this section first

% first, load a "freq" struct from any of our results
load('C:\Users\dzhu9728\Desktop\EEG\results_freq\9002_S1_EC.mat');

% change freq range & EEG channels
freq.freq = 4:25; % freq range in her data was 4~25Hz
freq.label = {'Fp1'; 'Fp2'; 'F7'; 'F3'; 'Fz'; 'F4'; 'F8';
    'FC5';'FC1';'FC2';'FC6';'T3';'C3';'Cz';'C4';'T4';'CP5';'CP1';
    'CP2';'CP6';'T5';'P3';'Pz';'P4';'T6';'O1';'O2'};

% cut matrix to correct size (22 freqs * 27 channels)
freq.powspctrm = freq.powspctrm(1:22, :);
freq.powspctrm = freq.powspctrm(:, 1:27);


%% Do some manual steps here, then run the rest

% <MANUAL STEP 1> paste the data into freq.powspctrm
%
% Flavia's EEG data in an Excel file:
% Z:\Lab_Manuscripts_Ethics_Grants\Manuscripts\Published\Published2018\Flavia_EEG_TCD\temp\send to Luke 19 vs 17 10052017


% <MANUAL STEP 2> specify the subject ID (plots for each subject will be saved into a separate folder)
subj_ID = '610';


% should run automatically from here onwards!

save_location = ['C:\Users\dzhu9728\Desktop\EEG\Flavia_results\' subj_ID '\'];
mkdir(save_location);

% transpose the matrix (FT plot rqs chan * freq)
freq.powspctrm = freq.powspctrm'; 

% plot power spectrum for all channels (overlay)
figure; hold on;
for chan = 1:length(freq.label)
    plot(freq.freq, freq.powspctrm(chan,:))
end
xlim([4 25]);
xlabel('Frequency (Hz)');
ylabel('Absolute power (uV^2)');
hold off;

export_fig(gcf, [save_location 'powspctrm_allchans.png']); % use this tool to save the figure exactly as shown on screen

% plot avg of all channels (log transformed)
figure; plot(freq.freq, mean(log(freq.powspctrm)));
xlim([4 25]);
xlabel('Frequency (Hz)');
ylabel('Power (log[uV^2]');

export_fig(gcf, [save_location 'powspctrm_avg.png']); % use this tool to save the figure exactly as shown on screen

% plot topography for each freq band
plot_TFR_topo(freq, lay, 'theta', [4 8], save_location)
plot_TFR_topo(freq, lay, 'alpha', [9 12], save_location)
plot_TFR_topo(freq, lay, 'beta', [13 25], save_location)
