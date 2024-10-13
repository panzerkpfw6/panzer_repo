clear all
close all
pwd
fname=['snap_fwd_6867.raw'];
fname=['snap_fwd_6803.raw'];
fname=['snapshot_10'];
fname=['snapshot_2000'];
% fname=['snapshot_5'];
dims=[2301,512,751];
dims=[2309,520,759];
ccnt=dims(1)*dims(2)*dims(3);

% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);

s = dir(fname);         
filesize = s.bytes 

domain_size=2301*512*751*4;
format longG
disp(domain_size)

data = reshape(A,dims);
min(data,[],'all')
max(data,[],'all')

figure
imagesc( squeeze(data(:,310,:)).' );
% imagesc( squeeze(data(310,:,:)).' );
% % imagesc( squeeze(data(:,:,310)).' );
val=1e-12;
clim([-val,val])
colorbar

ss=1;