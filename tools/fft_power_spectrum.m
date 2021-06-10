cfg = [];
cfg.output = 'pow';
cfg.pad    = 'nextpow2';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = 1;
cfg.keeptrials = 'no';
%cfg.channel = {1:160};
datapow = ft_freqanalysis(cfg, all_blocks); % should pass in epoched data here
                                            % (continuous data is usually too large 
                                            % for ft_freqanalysis to handle)

% plot the power spectrum
cfg         = [];
cfg.channel = {'AG081', 'AG085', 'AG089'}; % select channels here

figure;
ft_singleplotER(cfg, datapow);
grid minor;
xlabel('Frequency (Hz)');
title('Frequency Power Spectrum');
set(gca,'FontSize',10);
