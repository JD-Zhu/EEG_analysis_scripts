% Manually select which ICA components to reject:
% eye blinks, other large muscle artefacts (e.g. jaw clenching, hand movements)

function [data_clean] = ICA_reject_comps(data, comp, lay, output_path)

    %% Plot the ICA components, so we can identify which comps to remove

    % change the colourmap
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64, 'RdBu'))) % change the colormap
    
    nComps = length(comp.label); % how many ICA comps were created?
    if nComps > 49
        nComps = 49; % only plot the first 49 comps (fits on 1 page)
    end
    
    figure; % this line is necessary, or else ft_topoplotIC won't plot
    cfg = [];
    cfg.component = 1:nComps;   
    cfg.layout = lay;
    cfg.marker = 'off';
    cfg.comment = 'no';
    ft_topoplotIC (cfg, comp);

    set(gcf, 'Position', get(0, 'Screensize')); % make the figure full-screen
    export_fig(gcf, [output_path 'ICA_comps.png']); % save the ICA topoplot

    cfg          = [];
    %cfg.channel  = 1:10; % components to be plotted
    cfg.viewmode = 'component';
    cfg.blocksize = 60; % display 60-sec segments
    cfg.layout   = lay;
    ft_databrowser(cfg, comp);
    
    
    %% For quality check: plot the continuous data BEFORE rejecting comps
    cfg          = [];
    cfg.channel  = 'all'; % components to be plotted
    cfg.viewmode = 'vertical';
    cfg.blocksize = 120; % display 120-sec segments
    cfg.ylim     = [ -64  64 ];
    cfg.layout   = lay;
    ft_databrowser(cfg, data);
    
    drawnow;
    
    
    %% Remove components
    diary on;
    success = false;
    
    while ~success
        % Ask user to specify the components to be removed
        disp('Enter components in the form [1 2 3]'); drawnow;
        comps2remove = input('Which components would you like to remove?\n'); drawnow;

        % Remove these components
        cfg           = [];
        cfg.component = comps2remove; %these are the components to be removed
        data_clean    = ft_rejectcomponent(cfg, comp, data);


        % For quality check: plot the continuous data AFTER rejecting comps
        cfg          = [];
        cfg.channel  = 'all'; % components to be plotted
        cfg.viewmode = 'vertical';
        cfg.blocksize = 120; % display 120-sec segments
        cfg.ylim     = [ -64  64 ];
        cfg.layout   = lay;
        ft_databrowser(cfg, data_clean);
        drawnow;

        
        % Let user decide whether they want to retry with diff comps
        prompt = ['\nCompare the datasets with and without IC removal:\n' ...
            'If it is OK, press Y to continue.\nIf not ok, press N, ' ...
            'then you will be asked to select the components again.\n\n'];

        answer = input(prompt, 's'); % read keyboard input as a string
        if strcmp(answer, 'N')
            success = false;
        else % for any other key press, we treat it as 'Yes'
            success = true;
        end
    end
    
    diary off;

end