clear all
close all
pwd
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/salt3d_676x676x201_zyx.raw'];
% fname=['/Volumes/ssd1/SIMWAVE/simwave/models/salt3d_676x676x201_xyz.raw'];
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/marmousi_2301x512x751_xyz.raw'];
fname=['../../velocity_models/salt3d_676x676x201_xyz.raw'];

% dims=[2301,512,751];
dims=[676,676,201];
% dims=[201,676,676];
ccnt=dims(1)*dims(2)*dims(3);

% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);

format longG
s = dir(fname);         
filesize = s.bytes 
disp(strcat('theoretical file size= ',int2str(ccnt*4)))
disp(strcat('true file size= ',int2str(filesize)))

data = reshape(A,dims);
min(data,[],'all')
max(data,[],'all')

figure
imagesc( squeeze(data(:,310,:))' );
colorbar

ss=1;