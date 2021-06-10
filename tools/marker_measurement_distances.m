marker_files = dir('*.mrk');
mrkp         = zeros(length(marker_files),5,3);%preallocate: 5 marker coils, 3 coords (x y z)

% each cycle is one file
for i =1:length(marker_files)
    mrk         = ft_read_headshape(marker_files(i).name,'format','yokogawa_mrk');
    mrk         = ft_convert_units(mrk,'mm');
    mrkp(i,:,:) = mrk.fid.pos;
end

mrke = zeros(size(mrkp,2),1);%preallocate how many marker coils

% each cycle is one marker coil
for i = 1:size(mrkp,2)
    mrke(i) = mean(pdist(squeeze(mrkp(:,i,:)))); % mean movement for this marker coil
end

figure;bar(mrke)
ylabel('mean error (mm)')
xlabel('marker')

disp(['total mean error = ', num2str(mean(mrke)), 'mm'])

mean_coords=squeeze(mean(mrkp,1)); % mean coords for each coil (across all measurements)
