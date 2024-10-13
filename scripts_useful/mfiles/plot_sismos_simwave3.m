clear all
close all

fname=['../simwave_export_to_ecrc_servers/data/sismos_52.raw'];
fname=['../simwave_export_to_ecrc_servers/data/sismos_53.raw'];
fname=['../simwave_export_to_ecrc_servers/data/sismos_54.raw'];
fname=['../simwave_export_to_ecrc_servers/data/sismos_3333.raw'];
dims=[51,51,100];
fname=['/Volumes/ssd1/SIMWAVE/simwave/seismic/sismos_6865.txt'];
dims=[11730,100];
dims=[11730,112];

fname=['/Volumes/ssd1/snapshots/sismos_1301.raw'];
dims=[11730,100];
dims=[11730,112];
dims=[2601,1000]


ccnt=dims(1)*dims(2);

% some file size calculations
s = dir(fname);         
filesize = s.bytes 
domain_size=ccnt*4;
format longG
disp(strcat('theoretical file size= ',int2str(domain_size)))
disp(strcat('true file size= ',int2str(filesize)))


% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
% A = fscanf(fileID,'%f',ccnt);
fclose(fileID);
data = reshape(A,dims);

% plotting
figure
% imagesc( squeeze(data(:,:,11)).' );
imagesc(data)
val=1e-12;
% clim([-val,val])
colorbar

% figure
% imagesc( squeeze(data(:,25,:)).' );
% val=1e-12;
% % clim([-val,val])
% colorbar

min(data,[],'all')
max(data,[],'all')

ss=1