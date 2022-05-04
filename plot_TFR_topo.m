% Plots the topography of each frequency band &
% saves the figures into the specified location
%
% @freq             output struct from ft_freqanalysis or ft_freqgrandaverage
% @lay              channel layout
% @freq_band        string, e.g. 'theta'
% @freq_range       e.g. [4 8]
% @save_location    where to save the figures (optional).
%                   if you don't want to save, just specify as ''
% @varargin{1}      axis limits for t-values, e.g. [0 3]
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function plot_TFR_topo(freq, lay, freq_band, freq_range, save_location, varargin)

    cfg = [];
    cfg.xlim         = freq_range; % freq range to plot (NOTE this is diff from FT documentation, as we only have 2 dimensions in "freq", ie. not TFR)
    if length(varargin) > 0
        cfg.zlim     = varargin{1};
    end
    cfg.marker       = 'on';
    cfg.style        = 'straight'; % 'both';
    cfg.colorbar     = 'yes';
    cfg.comment      = 'no';
    %cfg.comment      = 'xlim'; % display date & freq range
    cfg.layout       = lay;
    cfg.dataname     = freq_band; % display this in the figure window title bar
    ft_topoplotTFR(cfg, freq); 

    if ~strcmp(save_location, '')
        %saveas(gcf, [save_location freq_band '.png']); % this fn does not maintain the aspect ratio, font size, etc
        export_fig(gcf, [save_location freq_band '.png']); % use this tool to save the figure exactly as shown on screen
    end
end