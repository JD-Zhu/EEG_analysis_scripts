ID = '677';
data_path = ['Z:\Analysis\Judy\EpisodicMigraine\data\migraineurs\']; % this directory should contain all the SubjectFolders

% find the raw EDF file
SubjectID = ['Subject_' ID];
folder = [data_path SubjectID '\'];
confile_name = '*\*.edf';
files = dir(fullfile(folder, confile_name));
rawfile = fullfile(files(1).folder, files(1).name);

% read in EC portion (first 290 sec)
EEG = pop_biosig(rawfile, 'blockrange',[0 290]); 
EEG.setname = ID;
EEG = eeg_checkset( EEG );

% read chanlocs & remove non-EEG channels
if length(EEG.chanlocs) == 33
    EEG = pop_chanedit(EEG, 'load',{'Z:\\Analysis\\Judy\\EEG_analysis_scripts\\chanlocs_XYZ_33chan_forEEGLAB.txt' 'filetype' 'sfp'});
elseif length(EEG.chanlocs) == 32
    EEG = pop_chanedit(EEG, 'load',{'Z:\\Analysis\\Judy\\EEG_analysis_scripts\\chanlocs_XYZ_32chan_forEEGLAB.txt' 'filetype' 'sfp'});
end
EEG = pop_select( EEG,'nochannel',{'EOG1' 'EOG2' 'REF' 'ECG1' 'ECG2' 'TRIG'});
EEG = eeg_checkset( EEG );

% plot channel spectra
figure('name',[SubjectID ' - spectopo()']); 
pop_spectopo(EEG, 1, [0  289998], 'EEG' , 'freq', [10 30 50], 'freqrange',[2 60],'electrodes','off');
EEG = eeg_checkset( EEG );


% Automatic channel rejection (Kurtosis method)
%EEG_Kurt = pop_eegfiltnew(EEG, 1,35,1650,0,[],1); % apply HPF & LPF first
%EEG_Kurt = pop_rejchan(EEG_Kurt, 'elec',[1:27] ,'threshold',5,'norm','on','measure','kurt');

% Automated artifact rejection with "Clean Rawdata" plugin
% https://eeglab.org/tutorials/06_RejectArtifacts/cleanrawdata.html
%EEG_ASR = pop_clean_rawdata(EEG, 'FlatlineCriterion','off','ChannelCriterion',0.8,'LineNoiseCriterion',4,'Highpass',[0.2 1.8] ,'BurstCriterion',20,'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
%pop_eegplot(EEG_ASR, 1, 1, 1);
