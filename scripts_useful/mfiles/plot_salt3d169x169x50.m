clear all
close all
pwd
fname='salt3d_169x169x50_xyz.raw';
dims=[169,169,50];
ccnt=dims(1)*dims(2)*dims(3);

% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);

% read file to the end of file
% fid = fopen(fname,'r');
% while ~feof(fid)
%     A2=fread(fid,ccnt,'float32');
% end
% fclose(fid);

% A2=load(fname);
% A3=textscan(fname,'%f');
% A4=textscan(fname,'%d','HeaderLines',9);

data = reshape(A,dims);
min(data,[],'all')
max(data,[],'all')

figure
imagesc( squeeze(data(:,:,40)).' );
colorbar

ss=1;