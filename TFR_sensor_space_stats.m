% MEG125_into_segments;
% 
% rois=[rfrontal_idx';
%         rtemporal_idx';
%         roccipital_idx';
%         lfrontal_idx';
%         ltemporal_idx';
%         loccipital_idx'];


% labels = {'theta','alpha','beta'};
% freqbands = {[4 7],[8 13],[12 26]};

% Delta wave ? (0.1 ? 3 Hz)
% Theta wave ? (4 ? 7 Hz)
% Alpha wave ? (8 ? 15 Hz)
% Mu wave ? (7.5 ? 12.5 Hz)
% SMR wave ? (12.5 ? 15.5 Hz) Sensorimotor rhythm e.g. Behavioral management of epileptic seizures following EEG biofeedback training of the sensorimotor rhythm
% Beta wave ? (16 ? 31 Hz)
% Gamma wave ? (32 ? 100 Hz)

%e.g., beta1: 13?18?Hz; beta2: 19?25?Hz; beta3: 26?30?Hz)

time_1 = -0.75;
time_2 = 0.0;
time_3 = 0.75;

freq_1 = 13;
freq_2 = 20;

% for d=1:size(rois,1)
%for d=2

orig      = cd;
folders   = dir;
data_pre  = [];
data_post = [];
for i=1:length(folders)
    if ~isempty(str2num(folders(i).name))
        cd([orig,'/',folders(i).name,'/MEG'])
        if exist('TFR.mat','file')
            load ('TFR.mat')
            
            time_1_idx = find(ismember(round(norm_TFR_broad.time*1000)/1000,time_1));
            time_2_idx = find(ismember(round(norm_TFR_broad.time*1000)/1000,time_2));
            time_3_idx = find(ismember(round(norm_TFR_broad.time*1000)/1000,time_3));
            
            N=norm_TFR_broad.freq;
            V                       = [N' repmat(freq_1,length(N),1) repmat(freq_2,length(N),1)];
            [minValue(1),freq_1_idx] = min(abs(V(:,1)-V(:,2)));
            [minValue(2),freq_2_idx] = min(abs(V(:,1)-V(:,3)));
            %closestValue            = V(closestIndex,1);
            
            eval(['TFR_',folders(i).name,'=norm_TFR_broad;'])
        % data_pre   = cat(3,data_pre,squeeze(mean(norm_TFR_broad.powspctrm(rois(d,:),:,time_1_idx:time_2_idx))));
          %  data_post  = cat(3,data_post,squeeze(mean(norm_TFR_broad.powspctrm(rois(d,:),:,time_2_idx:time_3_idx))));
%         %    
             data_pre  = cat(3,data_pre,squeeze(mean(norm_TFR_broad.powspctrm(:,freq_1_idx:freq_2_idx,time_1_idx:time_2_idx),2))); %chan X freq X time - mean over freq
             data_post = cat(3,data_post,squeeze(mean(norm_TFR_broad.powspctrm(:,freq_1_idx:freq_2_idx,time_2_idx:time_3_idx),2)));
            
%             data_pre  = cat(4,data_pre,(norm_TFR_broad.powspctrm(:,freq_1_idx:freq_2_idx,time_1_idx:time_2_idx))); %chan X freq X time - mean over freq
%             data_post = cat(4,data_post,(norm_TFR_broad.powspctrm(:,freq_1_idx:freq_2_idx,time_2_idx:time_3_idx)));
        else
        end
    else
    end
    cd(orig)
end

%data_pre=permute(data_pre,[4 1 3 2]); %function wants it as participants x channels x time x frequency
%data_post=permute(data_post,[4 1 3 2]);
% 
 data_pre=permute(data_pre,[3 1 2]); %function wants it as participants x  time x frequency
 data_post=permute(data_post,[3 1 2]); %function wants it as participants x  time x frequency


%V = randi(10,[5 1])
%N = randi(10,[5 1])


% A = repmat(N,[1 length(V)]);
% [minValue,closestIndex] = min(abs(A-V'));
% closestValue = N(closestIndex);



newlocs     = fieldtripchan2eeglab_paul(norm_TFR_broad.grad);
ept_tfce_nb = ept_ChN2(newlocs,1);

%[Input]
% Participant data should be organised into two factors
%     e.g. Dataset 1 is controls while Dataset 2 is patients
% Electrode locations file created using eeglab is required to
% % calculate channel neighbours
%
% Analysis Types
% i = independent sample T-Test
% d = dependent (paired) sample T-Test
% o = one-sample T-test
% c = Pearson's correlation (r)

%Participants x Channels x Samples
%Results = ept_TFCE(DataFile1, DataFile2, ElecFile, varargin)

%Results = ept_TFCE(data_post, data_pre, [], 'flag_ft',1,'plots',1,'type','d');


subjects             = importdata('/Users/mq20096022/Desktop/Paul_MEG/Child MEG data/subject_key.csv'); %importdata('myfile.txt', ' ', 1)
subjects.textdata(1) = subjects.textdata(2);

%// Get a list of everything in this folder
D = dir;

%// Grab just the directories and remove '.' and '..'
folders = {D([D.isdir]).name};
folders = folders(~ismember(folders, {'.', '..'}));

%// Get all the files
files = {D(~[D.isdir]).name};

%// Remove .DS_store files
files                                                                                = files(~strcmpi(files, '.DS_store'));
groups                                                                               = subjects.textdata;
groups(find(subjects.data==setdiff(subjects.data,str2num(cell2mat(folders')))))      = [];
groups(1)                                                                            = groups(2);
subjects_no                                                                          = subjects.data;
subjects_no(find(subjects.data==setdiff(subjects.data,str2num(cell2mat(folders'))))) = [];
[B,I]                                                                                = sort(subjects_no);
subject_assignment                                                                   = groups(I);

Results = ept_TFCE(data_post(ismember(subject_assignment,[{'a'}]),:,:), data_pre(ismember(subject_assignment,[{'a'}]),:,:), newlocs, 'flag_ft',0,'plots',1,'type','d');

% ept_TFCE_ANOVA({data_pre(ismember(subject_assignment,{'s'}),:,:), data_post(ismember(subject_assignment,{'s'}),:,:); data_pre(ismember(subject_assignment,{'c'}),:,:)...
%     ,data_post(ismember(subject_assignment,{'c'}),:,:); data_pre(ismember(subject_assignment,{'a'}),:,:), data_post(ismember(subject_assignment,{'a'}),:,:)}, newlocs)

% cfg = [];
% cfg.chan = 1:125; % the other channels are eye-movement channels from the Eye-ufresult toolbox
% %res = ept_TFCE(data_post(ismember(subject_assignment,{'s'}),:,:), data_post(ismember(subject_assignment,[{'c'}]),:,:),'nperm',500,'neighbours',ept_tfce_nb(cfg.chan,:),'flag_ft',0,'plots',1,'type','i');
% %[res,info] = ept_TFCE(struct('nperm',500,'neighbours',ept_tfce_nb(cfg.chan,:)),permute(data,[3 1 2]));
%res = ept_TFCE(data_post(ismember(subject_assignment,{'s'}),:,:,:), data_post(ismember(subject_assignment,[{'c'}]),:,:,:),newlocs,'nperm',500,'flag_ft',0,'plots',1,'type','i');
% 
% fgPlot = [];
% cfgPlot.pvalues = res.P_Values;
% cfgPlot.colormap = {{'div','RdYlBu'},{'div','RdBu'},'seq'};
% cfgPlot.topoalpha = 0.0005; % where to put the significance dots?
% cfgPlot.individualcolorscale = 'row'; % different rows have very different interpretation
% cfgPlot.time = [0 0.3]; % we will zoom in
%plot_topobutter(cat(3,mean(data(1:125,:,:),3),res.TFCE_Obs,res.P_Values),ufresult.times,ufresult.chanlocs(1:125),cfgPlot)

%end
