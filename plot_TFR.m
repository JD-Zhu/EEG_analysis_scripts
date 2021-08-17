% Plots the power spectrum & topography of each frequency band
%
% @freq             output struct from ft_freqanalysis or ft_freqgrandaverage
% @lay              channel layout
% @save_location    where to save the figures (optional) 
%                   if you don't want to save, just specify as ''
% @x_limits         xlim for the power spectrum plot
% @infraslow        whether to plot infraslow band (true or false)
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function plot_TFR(freq, lay, save_location, x_limits, infraslow)

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
    if infraslow
        plot_TFR_topo(freq, lay, 'infraslow', [0.03 0.06], save_location)
    end
    plot_TFR_topo(freq, lay, 'theta', [4 8], save_location)
    plot_TFR_topo(freq, lay, 'alpha', [9 12], save_location)
    plot_TFR_topo(freq, lay, 'beta', [13 25], save_location)
end
