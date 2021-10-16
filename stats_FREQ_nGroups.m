%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_FREQ_nGroups.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Statistical analysis of frequency results:
% Comparison across 3 groups (e.g. migraine phases / frequency of attacks)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PLEASE SPECIFY the folder for this statistical analysis
stats_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\migraine_phases\';
%stats_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\migraine_frequency\';

% PLEASE SPECIFY the subject groups for comparison
groups = {'GA_prodrome', 'GA_postdrome', 'GA_interictal'};
%groups = {'GA_lessThan1day', 'GA_1-2days', 'GA_moreThan3days'};
% note: you need to have these GA folders ready inside the stats_folder

% Also set up the figure legends (MAKE SURE the order here is same as in the "groups" variable above)
figure_legends = {'Prodrome', 'Postdrome', 'Interictal'};
%figure_legends = {'< 1 day', '1-2 days', '> 3 days'};

load('lay_NeuroPrax32.mat');


%% plot overall power (one line == one group)

x_limits = [2 30];

% plot the GA for each group
figure; hold on;
for g = 1:length(groups)
    load([stats_folder groups{g} '\GA_avg.mat']); % GA for this group
    plot(GA_freq.freq, mean(GA_freq.powspctrm));
end

xlim(x_limits);
xlabel('Frequency (Hz)');
ylabel('Absolute power (uV^2)');
legend(figure_legends);
hold off;

export_fig(gcf, [stats_folder 'overall_power_for_each_group.png']);


%% individual channel analysis: ANOVA at each channel for each freq (27 x 30 = 810 comparisons)
% (see Figure 2a in Flavia paper)

load([stats_folder groups{1} '\GA_individuals.mat']); % load one group as an example
N_chan = size(GA_freq_indi.powspctrm, 2); % number of channels
N_freq = size(GA_freq_indi.powspctrm, 3); % number of freqs

% initialise "chan x freq" matrix to store f-values
f_values = zeros(N_chan, N_freq);
p_values = zeros(N_chan, N_freq);

for i = 1:N_chan % loop through each channel
    for j = 1:N_freq % loop through each freq
        
        % collate data for anova (each subject should have a single value)
        % https://au.mathworks.com/help/stats/anova1.html
        % https://au.mathworks.com/matlabcentral/answers/159736-how-do-i-perform-unbalanced-anova
        data_for_anova = [];
        grouping_var = {};
        for g = 1:length(groups)
            load([stats_folder groups{g} '\GA_individuals.mat']);
            data = GA_freq_indi.powspctrm(:,i,j); % each subject has a single value 
                                                  % representing the power at this chan & this freq
            data_for_anova = vertcat(data_for_anova, data); % append the data from all subjects in this group
            % also append this group's label accordingly (same number of times as how many subjects were appended)
            for count = 1:length(data)
                grouping_var = [grouping_var; groups{g}]; 
            end
        end
        
        [p,tbl,stats] = anova1(data_for_anova, grouping_var, 'off');
        %t_values(i,j) = stats.tstat; % store the t-value into the "chan x freq" matrix
        f_values(i,j) = cell2mat(tbl(2,5));
        p_values(i,j) = cell2mat(tbl(2,6));
    end
end

save([stats_folder 'indi-chan-analysis\stats--chan_x_freq.mat'], 'f_values', 'p_values');

%% plots
% using a separate section here as the ANOVAs above take a while to run,
% so we just load the stats results prevly saved
load([stats_folder 'indi-chan-analysis\stats--chan_x_freq.mat']);

figure;
imagesc(f_values)
colorbar
ylabel('EEG channel');
xlabel('Frequency (Hz)');
export_fig(gcf, [stats_folder 'indi-chan-analysis\indi-chan-analysis_F-values.png']);

% not necessary to plot the p-values (it would just be the opposite of the f-values)
%{
figure;
imagesc(p_values, [0 0.05]) % only plot p-values up to 0.05
colorbar
ylabel('EEG channel');
xlabel('Frequency (Hz)');
export_fig(gcf, [stats_folder 'indi-chan-analysis\indi-chan-analysis_p-values.png']);
%}


%% ANOVA at each channel for each freq band (27 x 3 = 81 comparisons)
% (see Figure 2b in Flavia paper)

% need to specify each freq band manually for now
freq_range = 9:12;
freq_band = 'alpha';

load([stats_folder groups{1} '\GA_individuals.mat']); % load one group as an example
N_chan = size(GA_freq_indi.powspctrm, 2); % get the number of channels

% initialise an array to store the f-value for each channel
f_values = zeros(N_chan, 1);
p_values = zeros(N_chan, 1);

for i = 1:N_chan % loop through each channel
    % collate data for anova (each subject should have a single value - see above)
    data_for_anova = [];
    grouping_var = {};
    for g = 1:length(groups)
        load([stats_folder groups{g} '\GA_individuals.mat']);
        data = GA_freq_indi.powspctrm(:,i,freq_range); % extract the power values at this channel (only for the selected freq range)
        data = mean(data,3); % take the avg over that freq range (each subject now ends up with a single value)
        data_for_anova = vertcat(data_for_anova, data); % append the data from all subjects in this group
        % also append this group's label accordingly (same number of times as how many subjects were appended)
        for count = 1:length(data)
            grouping_var = [grouping_var; groups{g}]; 
        end
    end

    [p,tbl,stats] = anova1(data_for_anova, grouping_var, 'off');
    f_values(i) = cell2mat(tbl(2,5));
    p_values(i) = cell2mat(tbl(2,6));
end


% create a dummy var for plotting
freq = GA_freq;
freq.powspctrm = f_values;
freq.freq = 0; % we no longer have a frequency dimension, just fill with a dummy value

% plot topography based on the f-values
plot_TFR_topo(freq, lay, freq_band, [], [stats_folder 'indi-chan-analysis\fvalues_']);

