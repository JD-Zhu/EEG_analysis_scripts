%
% Author: Judy Zhu (github.com/JD-Zhu)
%

% find all subject folders in the MRI database
MRI_database = 'D:\\Judy\\MRI_database\\';
SubjectFolders = listFolders(MRI_database);

nii_gz_filename = '\\MEG\\anatomy\\T1w_acpc_dc_restore.nii.gz'; % the zip file to extract

for i = 1:length(SubjectFolders)
    gunzip([MRI_database SubjectFolders{i} nii_gz_filename]);
end
