clear all
close all

root='/Volumes/ssd1/snapshots';
fname=['snapshot_2000'];
fname=['SB_2nd.raw'];
fname=['SB_2nd_2017.raw'];
fname=['SB_2nd_2517.raw'];

% % fname2=['SB_1st.raw'];
% fname2=['TB_2nd.raw'];
% fname2=['TB_1st.raw'];
fname2=['SB_2nd_abc.raw'];
% fname2=['SB_2nd_abc_1017.raw'];
% % fname2=['SB_2nd_abc_2017.raw'];
fname2=['SB_2nd_abc_2517.raw'];
% fname2=['SB_2nd_abc_2517_2.raw'];
% fname2=['SB_2nd_abc_4017.raw'];
% fname2=['SB_2nd_abc_217.raw'];
fname=fullfile(root,fname);
fname2=fullfile(root,fname2);

dims=[2301,512,751];
% add 8 points to each dimension
dims=[2309,520,759];
dims=[520,520,520];
ccnt=dims(1)*dims(2)*dims(3);

% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data = reshape(A,dims);

fileID = fopen(fname2,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data2 = reshape(A,dims);

ix=10;
% plotting
figure
imagesc( squeeze(data(:,:,ix)).' );
val=1e-12;
% clim([-val,val])
colorbar

figure
imagesc( squeeze(data2(:,:,ix)).' );
colorbar

figure
imagesc( squeeze(data2(:,:,260)).' );
colorbar

figure
imagesc( squeeze(data(:,:,ix)-data2(:,:,ix)).' );
colorbar

min(data,[],'all')
max(data,[],'all')

min(data2,[],'all')
max(data2,[],'all')

diff=data-data2;
min(diff,[],'all')
max(diff,[],'all')

% some file size calculations
s = dir(fname);         
filesize = s.bytes 
domain_size=2301*512*751*4;
format longG
disp(domain_size)
pp=1


%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% limits = [NaN NaN NaN NaN NaN 10]
% [x, y, z, D] = subvolume(data2, limits);           
% [fo,vo] = isosurface(x,y,z,D,5);               
% [fe,ve,ce] = isocaps(x,y,z,D,5);               
% figure
% p1 = patch('Faces', fo, 'Vertices', vo);       
% p1.FaceColor = 'red'
% ss=1
%%%%%%%%%%%%%%%%%%%%%%%%