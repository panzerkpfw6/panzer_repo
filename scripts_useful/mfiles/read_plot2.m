% [SIMWAVE MSG]:... simwave information:
% [SIMWAVE MSG]:... velocity size       = 2301 x 512 x 748
% [SIMWAVE MSG]:... velocity min/max    = 1500.000000 - 5500.000000
% [SIMWAVE MSG]:... compute domain size = 2301 x 512 x 748 (865.663605 MB)
% [SIMWAVE MSG]:... imaging domain size = 2301 x 512 x 748 (840.404297 MB)
% [SIMWAVE MSG]:... cdp,line,depth      = 10 x 10 x 10

clear all
close all
% models/marmousi_2301x512x751_xyz.raw
pwd
fname='./test_wavefield2/snap_fwd_6867.raw'
dims=[1,2301,512,748];

ccnt=dims(1)*dims(2)*dims(3);
ccnt=dims(2)*dims(3)*dims(4);
% ccnt=dims(1)*dims(2)*dims(3)*dims(4);
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
length(A)
length(A)/2301/748
% data = reshape(A,dims);
% data = reshape(A,[2301,512,248,10]);
% data = reshape(A,[dims(1),dims(2),dims(3),dims(4)]);
data = reshape(A,[dims(2),dims(3),dims(4)]);
min(data,[],'all')
max(data,[],'all')

figure
% imagesc( squeeze(data(22,:,:)).' );
imagesc( squeeze(data(:,22,:,3)).' );
colorbar

ss=1;