% Plots ERF for each condition. 
% Generates and plots GFP.

% Can use this to regen all plots from saved erf results:
%load([ResultsFolder SubjectID '_erf.mat']); % select which subject to load
%load([ResultsFolder 'lay.mat']); 

% If you want to compare before & after artefact removal (erf vs. erf_clean),
% set plot_uncleaned to 'true' & supply the 'erf' input arg;
% if only plotting cleaned erf & GFP, supply the 'erf' arg as []
%
% @param erf:           uncleaned erf (ie. w/o PCA, w/o visual rejection)
% @param erf_clean:     normal erf
% @param erf_combined:  avg erf across all conds
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
function plot_ERF (erf, erf_clean, erf_combined, lay, plot_uncleaned, plot_multiplot, plot_combined)

    % run the #define section
    %global conds_cue; global conds_target; 
    global eventnames_real; global colours_and_lineTypes; 
    global colours; global lineTypes;
    global PLOT_XLIM; global ERF_BASELINE;
    common();
    

    % required configurations b4 any calls to ft_multiplotER()
    cfg              = [];
    cfg.showlabels   = 'yes';
    cfg.fontsize     = 6;
    cfg.layout       = lay;
    cfg.baseline     = ERF_BASELINE; % makes no diff if we've already done baseline correction earlier
    cfg.baselinetype = 'absolute';
    cfg.graphcolor   = cell2mat(colours); 
    cfg.linestyle    = lineTypes;
    

    % all pairwise comparisons btwn raw data (erf) & after artefact removal (erf_clean)
    % just to see how good the artefact removal is (atm we reject components 1:5, we can adjust this for more/less removal)
    if (plot_uncleaned)
        for j = 1:length(eventnames_real)
            figure;
            ft_multiplotER(cfg, erf.(eventnames_real{j}), erf_clean.(eventnames_real{j}));
        end
    end


    % to compare the 4 conds
    % (requires the list of 'cfg' assignments above, if running in console)
    if (plot_multiplot)
        %figure('Name','ft_multiplotER: erf_clean.cuechstay, erf_clean.cuechswitch, erf_clean.cueenstay, erf_clean.cueenswitch');
        figure('Name','ft_multiplotER: erf_clean');
        cfg.xlim = PLOT_XLIM;

        % convert struct to cell array, then you can feed it in as 'varargin'
        cellarray = struct2cell(erf_clean);
        ft_multiplotER(cfg, cellarray{:});
        % this replaces the lengthy way of manually writing out the list of conds (below)
        %{
        if length(eventnames_real) == 9 % plot 9 conds (ie. collapsed across langs)
            ft_multiplotER(cfg, erf_clean.NatStay, erf_clean.NatSwitch, erf_clean.NatSingle, ... 
                erf_clean.ArtStay, erf_clean.ArtSwitch, erf_clean.ArtSingle, ...
                erf_clean.BiStay, erf_clean.BiSwitch, erf_clean.BiSingle);
        elseif length(eventnames_real) == 18 % plot 18 conds
            cfg.graphcolor   = colours; % temp solution of using default colours, coz the custom colours aren't working for 18 conds yet
            ft_multiplotER(cfg, erf_clean.NatStayC, erf_clean.NatStayE, erf_clean.NatSwitchC, erf_clean.NatSwitchE, erf_clean.NatSingleC, erf_clean.NatSingleE, ... 
                erf_clean.ArtStayC, erf_clean.ArtStayE, erf_clean.ArtSwitchC, erf_clean.ArtSwitchE, erf_clean.ArtSingleC, erf_clean.ArtSingleE, ...
                erf_clean.BiStayC, erf_clean.BiStayE, erf_clean.BiSwitchC, erf_clean.BiSwitchE, erf_clean.BiSingleC, erf_clean.BiSingleE);
        else % for partial exps (only for testing purpose)
            fprintf('plot_ERF (line 55): Warning: check the number of conds in eventnames_real. The array length should be either 9 or 18.\n');
            % !modify the cond list below manually to match your list!
            ft_multiplotER(cfg, erf_clean.NatStay, erf_clean.NatSwitch, erf_clean.NatSingle, ... 
            erf_clean.ArtSingle, erf_clean.BiSingle);
        end
        %}
        
        % specify the legends manually (otherwise it will display incorrectly)
        lines = findall(gcf, 'Type','line');
        lines = lines([29 26 23 20 17 14 11 8 5]); % grab the correct lines (this is complicated because many lines are plotted in ft_multiplot)
        legend(lines, eventnames_real);


        %{
        figure('Name','ft_multiplotER: erf_clean.targetchstay, erf_clean.targetchswitch, erf_clean.targetenstay, erf_clean.targetenswitch');
        cfg.xlim = [-0.1 0.75];
        ft_multiplotER(cfg, erf_clean.targetchstay, erf_clean.targetchswitch, erf_clean.targetenstay, erf_clean.targetenswitch);
        legend(eventnames_real(conds_target));
        %}
    end


    %% Calc global averages across all sensors (GFP = global field power)
    cfg        = [];
    cfg.method = 'power';
    %cfg.channel = alldata.label([81:89]); % trigger spike is present in all of these channels: 81-89
    for j = 1:length(eventnames_real)
        erf_clean_GFP.(eventnames_real{j}) = ft_globalmeanfield(cfg, erf_clean.(eventnames_real{j}));
    end

    % plot GFP for cue-locked 
    figure('Name','GFP'); hold on
    for j = 1:length(eventnames_real)
    %for j = conds_cue
        %{
        % epoch was [-1 1], we only want to plot [-0.2 0.75]
        % so here we customise the time interval to plot
        start_time = find(erf_clean_GFP.(eventnames_real{j}).time >= -0.2);
        start_time = start_time(1);
        end_time = find(erf_clean_GFP.(eventnames_real{j}).time <= 0.75);
        end_time = end_time(end);
        %}
        if colours_and_lineTypes % use a combination of colours and line types to distinguish conds
            plot(erf_clean_GFP.(eventnames_real{j}).time, erf_clean_GFP.(eventnames_real{j}).avg, 'color', colours{j}, 'LineStyle', lineTypes{j});
        else % just use colours
            plot(erf_clean_GFP.(eventnames_real{j}).time, erf_clean_GFP.(eventnames_real{j}).avg, 'color', colours(j,:));
        end
        
        xlim(PLOT_XLIM); % epoch was [-1 1], we only want to plot [-0.1 0.75]
    end
    legend(eventnames_real);

    %{
    % plot GFP for target-locked 
    figure('Name','GFP_target'); hold on
    for j = conds_target
        plot(erf_clean_GFP.(eventnames_real{j}).time, erf_clean_GFP.(eventnames_real{j}).avg);
        xlim([-0.1 0.75]); % epoch was [-1 1], we only want to plot [-0.1 0.75]
    end
    legend(eventnames_real(conds_target));
    %}
    
    
    %% Sanity check for AEF pilot:
    % plot the indi-channel ERFs avg'd over all conds
    % we can check if sound offset (at 400ms) generates any brain response,
    % esp. in the auditory cortex!!
    if (plot_combined)
        cfg              = [];
        cfg.showlabels   = 'yes';
        cfg.fontsize     = 6;
        cfg.layout       = lay;
        cfg.baseline     = ERF_BASELINE; % makes no diff if we've already done baseline correction earlier
        cfg.baselinetype = 'absolute';
        cfg.xlim = PLOT_XLIM;

        figure('Name','multiplot: avg over conds'); 
        ft_multiplotER(cfg, erf_combined);


        % Plot the GFP avg'd over all conds
        cfg        = [];
        cfg.method = 'power';
        GFP_combined = ft_globalmeanfield(cfg, erf_combined);

        figure('Name','GFP: avg over conds');
        plot(GFP_combined.time, GFP_combined.avg);
        xlim(PLOT_XLIM);
    end
    
end
