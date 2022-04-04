% Plots the power spectrum & topography of each frequency band
%
% @freq             output struct from ft_freqanalysis or ft_freqgrandaverage
% @lay              channel layout
% @save_location    where to save the figures (optional) 
%                   if you don't want to save, just specify as ''
% @x_limits         xlim for the power spectrum plot
% @freq_bands       which freq bands to plot (e.g. theta, alpha, beta)
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function plot_TFR(freq, lay, save_location, x_limits, freq_bands)

    % unsquaring (i.e. back to uV unit)
    %bu = freq.powspctrm;
    %freq.powspctrm = sqrt(freq.powspctrm);

    % plot power spectrum of a particular channel
    %{
    figure;
    plot(freq.freq, freq.powspctrm(12,:))
    xlabel('Frequency (Hz)');
    ylabel('absolute power (uV^2)');
    xlim([1 30]);
    %}

    % plot power spectrum for all channels (overlay)
    figure; hold on;
    for chan = 1:length(freq.label)
        plot(freq.freq, freq.powspctrm(chan,:))
    end
    xlim(x_limits);
    xlabel('Frequency (Hz)');
    ylabel('Absolute power (uV^2)');
    hold off;

    export_fig(gcf, [save_location 'powspctrm_allchans.png']); % use this tool to save the figure exactly as shown on screen

    % plot avg of all channels
    figure; plot(freq.freq, mean(freq.powspctrm));
    xlim(x_limits); %[0.01 30]
    xlabel('Frequency (Hz)');
    ylabel('Absolute power (uV^2)');
    
    export_fig(gcf, [save_location 'powspctrm_avg.png']); % use this tool to save the figure exactly as shown on screen

    % plot avg of all channels (log transformed)
    figure; plot(freq.freq, mean(log(freq.powspctrm)));
    xlim(x_limits); %[0.01 30]
    xlabel('Frequency (Hz)');
    ylabel('Power (log[uV^2]');

    export_fig(gcf, [save_location 'powspctrm_avg_log.png']); % use this tool to save the figure exactly as shown on screen

    %{
    figure; 
    cfg = [];
    cfg.layout = lay;
    ft_multiplotER(cfg, freq);
    %}

    % topoplot for each freq band
    for band = 1:length(freq_bands)
        freq_band = cell2mat(freq_bands{band, 1}); % first field is the freq band name
        freq_range = freq_bands{band, 2}; % second field is the freq range in Hz
        plot_TFR_topo(freq, lay, freq_band, [freq_range(1) freq_range(end)], save_location);
    end
end
