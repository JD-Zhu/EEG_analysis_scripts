% You can use this script to plot the source activity timecourse in a single subject.
% "Source activity" can be for an ROI (i.e. virtual sensor), or for a single vertex,
% just load the desired source timecourse into the VE variable
%
% Author: Judy Zhu (github.com/JD-Zhu)
%

% SELECT which ROI/vertex to plot
VE = ROI_activity.V1; % load the ROI output file (e.g. A01-XC-3489_ROI.mat) first
%VE = singleVertex50; % load the saved timecourses (e.g. singleVertex_timecourse_in_A01_LIFG.mat) first

% SELECT which conds to plot:
conds = 1:9; % plot all conds
%conds = [1 3 4 6]; 


eventnames_real = fieldnames(VE);

figure; hold on; 
for j = conds
   plot(VE.(eventnames_real{j}).time, VE.(eventnames_real{j}).avg);
   xlim([-0.8 1]);
end
legend(eventnames_real(conds));
