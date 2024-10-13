clear all
close all

root='/Volumes/ssd1/snapshots';
fname= ['SB-final_snap_504.raw'];
fname2=['TB-final_snap_504.raw'];
fname3=['diff_TB-final_snap_504.raw'];
fname=fullfile(root,fname);
fname2=fullfile(root,fname2);
fname3=fullfile(root,fname3);

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

fileID = fopen(fname3,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data3= reshape(A,dims);

data=permute(data, [3 2 1]);
data2=permute(data2, [3 2 1]);
data3=permute(data3, [3 2 1]);

close all
ix=500;
ix=5;
ix=260;
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
imagesc( squeeze(data3(:,:,ix)).' );
colorbar

ss=1