%
% Author: Judy Zhu (github.com/JD-Zhu)
%

% run the #define section to obtain values for global vars
global ResultsFolder; 
global colours; global lineTypes; global PLOT_XLIM;
global PLOT_SHADE; % for plotting shaded boundary on each time course
common();

% SELECT which set of single-subject ERFs to use
run_name = 'TSPCA10000_3'; % this should be a folder name inside the "Results_ERF" folder
ResultsFolder_thisrun = [ResultsFolder run_name '\\']; % ERF results for all subjects


% load the results
load([ResultsFolder_thisrun 'GA_erf.mat']);
if ~strcmp(PLOT_SHADE, 'no') % load indi GA only if plotting shaded boundary (otherwise don't load -> saves memory)
    load([ResultsFolder_thisrun 'GA_individuals.mat']);
end
load([ResultsFolder_thisrun 'stats.mat']);


%% Collapse GA into 2 conds:
% below is only for target-locked, can do the same for cue-locked
%{
% (Eng vs Chn)
en = GA_erf.targetenstay;
en.avg = (GA_erf.targetenstay.avg + GA_erf.targetenswitch.avg) / 2;
ch = GA_erf.targetchstay;
ch.avg = (GA_erf.targetchstay.avg + GA_erf.targetchswitch.avg) / 2;

% (stay vs switch)
st = GA_erf.targetchstay;
st.avg = (GA_erf.targetchstay.avg + GA_erf.targetenstay.avg) / 2;
sw = GA_erf.targetchswitch;
sw.avg = (GA_erf.targetchswitch.avg + GA_erf.targetenswitch.avg) / 2;

% baseline correction
cfg = [];
cfg.baseline = [-0.1 0];
en = ft_timelockbaseline(cfg, en); 
ch = ft_timelockbaseline(cfg, ch);
st = ft_timelockbaseline(cfg, st); 
sw = ft_timelockbaseline(cfg, sw);
%}


%% Effect 1: target_lang, 360-515ms, p = 0.02
stat = target_lang;
start_time = 0.360;
end_time = 0.515;

% produce a list of all channels that were sig 
% at one or more time points during the effect interval
sig_channels = [];

% each cycle checks one channel
for i = 1:size(stat.mask, 1)
    if find(stat.mask(i,:)) % this channel was sig at some point
        sig_channels = [sig_channels i]; % so add it to the list
    end
end


% collapse into 2 conditions, based on the cross-subject grand avg
% (these will be the 2 lines in the plot, each line shows the time course for 1 collapsed cond)
en = GA_erf.targetenstay;
en.avg = (GA_erf.targetenstay.avg + GA_erf.targetenswitch.avg) / 2;
ch = GA_erf.targetchstay;
ch.avg = (GA_erf.targetchstay.avg + GA_erf.targetchswitch.avg) / 2;

% baseline correction
cfg = [];
cfg.baseline = [-0.1 0];
en = ft_timelockbaseline(cfg, en); 
ch = ft_timelockbaseline(cfg, ch);


% if plotting shaded boundary, need to do the following:
% collapse into 2 conditions, based on each individual subject's ERF time course
% (these will be used to calc the shaded boundary around each line in the plot)
if ~strcmp(PLOT_SHADE, 'no')
    en_indi = GA_indi.targetenstay;
    en_indi.individual = (GA_indi.targetenstay.individual + GA_indi.targetenswitch.individual) / 2;
    ch_indi = GA_indi.targetchstay;
    ch_indi.individual = (GA_indi.targetchstay.individual + GA_indi.targetchswitch.individual) / 2;

    % for each individual subject, average over all sig channels
    en_indi.individual = nanmean(en_indi.individual(:,sig_channels,:), 2); % take the mean on the 2nd dimension (i.e. channel)
    ch_indi.individual = nanmean(ch_indi.individual(:,sig_channels,:), 2); 
end
    

% plot
figure('Name', 'Average ERF of significant channels: target_lang_360-515ms'); hold on;

if strcmp(PLOT_SHADE, 'no') % do not plot shaded boundary, only plot the lines  
    cfg        = [];
    cfg.channel = sig_channels; % only include the sig channels
    cfg.linewidth = 3;
    ft_singleplotER(cfg, en, ch); % this autoly avg over the selected channels
                                  % alternatively, you can calc the mean manually (see under "else")
else
    margin_en = calc_margin(en_indi.individual, PLOT_SHADE);
    margin_ch = calc_margin(ch_indi.individual, PLOT_SHADE);

    % for the GA across subjects, average over all sig channels
    % (this should give you the same lines that ft_singleplotER autoly produces)
    en_avg = nanmean(en.avg(sig_channels,:), 1); % take the mean on the 1st dimension (i.e. channel)    
    ch_avg = nanmean(ch.avg(sig_channels,:), 1); 

    % plot time courses with shaded boundary
    boundedline(en.time, en_avg, margin_en(:), 'alpha', 'transparency',0.15, colours(1));
    boundedline(ch.time, ch_avg, margin_ch(:), 'alpha', 'transparency',0.15, colours(2));        
end

xlim([-0.1 0.55]);
ylim([-7e-15 7e-15]);
xlabel('Seconds');
ylabel('Tesla');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded

% specify the legend manually (otherwise it will include
% each shaded patch as an item too). For some reason,
% the order of the lines are reversed when you grab them
lines = findall(gcf, 'Type','line');
legend(flip(lines), {'English (L2)', 'Mandarin (L1)'}, 'Location','northwest', 'FontSize',30);

set(lines, 'Linewidth',3); % line thickness
hold off;


%% Effect 2: target_ttype, 355-465ms, p = 0.07 (marginal)
% the code below only plots single lines atm, no shaded boundary

stat = target_ttype;
start_time = 0.355;
end_time = 0.465;

% Extra step for marginal effects:
% Because the "stats" output we saved had alpha=0.05, so marginal effect 
% is not showing up in the .mask field. Here we recreate the .mask field 
% to show marginal effect.
% This approach for creating the .mask field has been verified to be correct.
idx = find(stat.prob < 0.1); % find indices where p < 0.1
stat.mask(idx) = 1; % fill these positions with 1s


% produce a list of all channels that were sig 
% at one or more time points during the effect interval
sig_channels = [];

% each cycle checks one channel
for i = 1:size(stat.mask, 1)
    if find(stat.mask(i,:)) % this channel was sig at some point
        sig_channels = [sig_channels i]; % so add it to the list
    end
end

% plot the ERF (averaged over all sig channels)
figure('Name', 'Average ERF of significant channels: target_ttype_355-465ms'); hold on;
cfg        = [];
cfg.channel = sig_channels; % only include the sig channels
cfg.linewidth = 3;
ft_singleplotER(cfg, st, sw);

%plot(st.time, st.avg, 'LineWidth',3);
%plot(sw.time, sw.avg, 'LineWidth',3);
xlim([-0.1 0.75]);
xlabel('Seconds');
ylabel('Tesla');
set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
legend({'Stay', 'Switch'}, 'Location','northwest', 'FontSize',30);
box on; % draw a border around the figure

% create shaded region indicating effect duration
ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
y = [ylow ylow yhigh yhigh];
patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
ylim(ylimits); % ensure ylim doesn't get expanded
hold off;
