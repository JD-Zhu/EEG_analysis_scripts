% Script to apply smoothing on all single-subject ROI timecourses, 
% to make them less "spiky" (could help to detect effects more easily)
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% 
%% setup
global ResultsFolder_ROI; % all subjects' ROI data are stored here
common();

% SPECIFY location containing the original (ie. not yet smoothed) timecourses
original_folder = [ResultsFolder_ROI 'TSPCA10000_3_freeori'];
% SPECIFY smoothing method/parameters: 'basic-smooth', 'median-3', 'median-4', 'SG-example', 'SG-paper', 'SG-paper-adjusted', 'SG-paper-adjusted-less'  
smooth_method = 'median-3';


%% start
% create a folder to store the smoothed timecourses
smoothed_folder = [original_folder '_' smooth_method '\\'];
mkdir(smoothed_folder);
original_folder = [original_folder '\\'];

% find all .mat files in original_folder
files = dir([original_folder '*_ROI.mat']);

% each cycle reads in one '.mat' file (ie. one subject's ROI timecourses)
for i = 1:length(files)
    filename = [original_folder files(i).name];
    load(filename);
    
    % loop thru each ROI
    ROIs_label = fieldnames(ROI_activity);
    for k = 1:length(ROIs_label)
        ROI_name = ROIs_label{k};
        
        % loop thru each cond
        conds_label = fieldnames(ROI_activity.(ROI_name));
        for j = 1:length(conds_label)
            cond_name = conds_label{j};
            
            % read the original (ie. not yet smoothed) timecourse
            timecourse = ROI_activity.(ROI_name).(cond_name).avg;
            
            % smooth it (using one of the options below)
            % https://www.ncbi.nlm.nih.gov/pubmed/32024267
            switch smooth_method
                case 'basic-smooth' % basic smooth fn
                    timecourse_smooth = smooth(timecourse); 
                    timecourse_smooth = timecourse_smooth'; % the basic smooth fn outputs a transposed vector, so we flip it back
                case 'median-3' % median filter with window size of 3
                    timecourse_smooth = medfilt1(timecourse, 3); % adjusted window sizes based on sampling rate (theirs is 512Hz, mine is downsampled to 200Hz)
                case 'median-4'
                    timecourse_smooth = medfilt1(timecourse, 4);
                case 'SG-example' % https://www.mathworks.com/help/signal/ref/sgolayfilt.html
                    timecourse_smooth = sgolayfilt(timecourse, 3, 11);
                case 'SG-paper' % see ref paper above
                    timecourse_smooth = sgolayfilt(timecourse, 4, 27);
                case 'SG-paper-adjusted' % adjusted window size based on sampling rate (theirs is 512Hz, mine is downsampled to 200Hz)
                    timecourse_smooth = sgolayfilt(timecourse, 4, 11); 
                case 'SG-paper-adjusted-less' % just trying a bit larger window, hopefully combines 2 peaks into 1
                    timecourse_smooth = sgolayfilt(timecourse, 4, 13); 
            end
            
            % replace the original
            ROI_activity.(ROI_name).(cond_name).avg = timecourse_smooth;
            
            % plot for sanity check
            %{
            figure;
            plot(ROI_activity.(ROI_name).(cond_name).time, timecourse);
            hold on;
            plot(ROI_activity.(ROI_name).(cond_name).time, timecourse_smooth);
            hold off;
            %}
        end
    end
    
    % save the smoothed timecourses to new location
    save([smoothed_folder files(i).name], 'ROI_activity');
end
