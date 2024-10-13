clear all
close all

dims=[2301,512,751];
% add 8 points to each dimension
dims=[2309,520,759];
fname=['snapshot_2000'];

dims=[520,520,520];
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/salt3d_512x512x512_zyx_pad.raw'];

dims=[520,520,520];
% fname=['/Volumes/ssd1/SIMWAVE/simwave/models/salt3d_512x512x512_zyx_pad.raw'];
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/salt3d_512x512x512_zyx_pad.raw'];

dims=[684,684,209];
dims=[676,676,201];
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/salt3d_676x676x201_xyz.raw'];

% dims=[520,520,520];
% fname=['/Volumes/ssd1/SIMWAVE/simwave/models/two_layer_1500_4500_zyx_pad.raw'];

% dims=[759,520,2309];
% dims=[751,512,2301];
% fname=['/Volumes/ssd1/SIMWAVE/marmousi_751x512x2301_zyx.raw'];

% dims=[512,512,512];
% fname=['/Volumes/ssd1/SIMWAVE/simwave/models/benchmark_512_3d.raw'];

dims=[2048,2048,512];
fname=['benchmark_2048_3d.raw'];

dims=[512,512,512];
fname=['benchmark_512_3d.raw'];


% some file size calculations
s = dir(fname);         
filesize = s.bytes ;
domain_size=2301*512*751*4;
domain_size=dims(1)*dims(2)*dims(3)*4;
format longG
disp(strcat('true size=',int2str(filesize)))
disp(strcat('calculated domain size=',int2str(domain_size)))

ccnt=dims(1)*dims(2)*dims(3);
% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
% A2=fread(fileID,[],'float32');
fclose(fileID);
data = reshape(A,dims);
% data = reshape(A,[520,520,520]);

% plotting
figure
imagesc( squeeze(data(:,100,:)).' );
val=1e-12;
% clim([-val,val])
colorbar

figure
imagesc( squeeze(data(:,:,170)).' );
colorbar
ss=1


