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
stats_folder = 'Z:\Analysis\Judy\EpisodicMigraine\stats\migraine_frequency\';

% PLEASE SPECIFY the subject groups for comparison
groups = {'GA_lessThan1day', 'GA_1-2days', 'GA_moreThan3days'};
% you need to have these GA folders ready inside the stats_folder

load('lay_NeuroPrax32.mat');


%% plot overall power (one line == one group)

% set up for the figure
legends = {'< 1 day', '1-2 days', '> 3 days'}; % MAKE SURE the order here is same as in the "groups" variable above
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
legend(legends);
hold off;

export_fig(gcf, [stats_folder 'overall_power_for_each_group.png']);


%% individual channel analysis: 25vs12 t-test at each freq for each channel
% (Figure 2a in Flavia paper)

load([stats_folder groups{1} '\GA_individuals.mat']); % load one group as an example
N_chan = size(GA_freq_indi.powspctrm, 2); % number of channels
N_freq = size(GA_freq_indi.powspctrm, 3); % number of freqs

% initialise "chan x freq" matrix to store f-values
f_values = zeros(N_chan, N_freq);
p_values = zeros(N_chan, N_freq);

for i = 1:N_chan % loop through each channel
    for j = 1:N_freq % loop through each freq
        
        % collate data for anova (each subject should have a single value 
        % representing the power at this chan & this freq)
        % https://au.mathworks.com/help/stats/anova1.html
        % https://au.mathworks.com/matlabcentral/answers/159736-how-do-i-perform-unbalanced-anova
        data_for_anova = [];
        grouping_var = {};
        for g = 1:length(groups)
            load([stats_folder groups{g} '\GA_individuals.mat']);
            data = GA_freq_indi.powspctrm(:,i,j);
            data_for_anova = vertcat(data_for_anova, data); % append the data for subjects belonging to this group
            % also append this group's label accordingly (same number of times as how many subjects were appended)
            for count = 1:length(data)
                grouping_var = [grouping_var; groups{g}]; 
            end
        end
        
        [p,tbl,stats] = anova1(data_for_anova, grouping_var, 'off');
        %t_values(i,j) = stats.tstat; % store the f-value into the "chan x freq" matrix
        f_values(i,j) = cell2mat(tbl(2,5));
        p_values(i,j) = cell2mat(tbl(2,6));
    end
end

save([stats_folder 'indi-chan-analysis\stats--chan_x_freq.mat'], 'f_values', 'p_values');

%%
figure;
imagesc(f_values)
colorbar
ylabel('EEG channel');
xlabel('Frequency (Hz)');
export_fig(gcf, [stats_folder 'indi-chan-analysis\indi-chan-analysis_F-values.png']);

figure;
imagesc(p_values, [0 0.05]) % only plot p-values up to 0.05
colorbar
ylabel('EEG channel');
xlabel('Frequency (Hz)');
export_fig(gcf, [stats_folder 'indi-chan-analysis\indi-chan-analysis_p-values.png']);

