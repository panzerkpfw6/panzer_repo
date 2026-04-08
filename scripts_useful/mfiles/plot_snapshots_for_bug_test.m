clear all
close all

root='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/PhD_proposal/work_with_david/Exawave_3_handover/stencil/scripts/python/bug_test';
fname= ['SB_2nd.raw'];
fname2=['TB_2nd.raw'];

fname=fullfile(root,fname);
fname2=fullfile(root,fname2);

% add 8 points to each dimension
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

data=permute(data, [3 2 1]);
data2=permute(data2, [3 2 1]);

close all

ix=260;
% plotting
figure
imagesc( squeeze(data(:,:,ix)).'  );
% val=1.5;clim([-val,val]);
title('SB 2nd');
colorbar

figure
imagesc( squeeze(data2(:,:,ix)).' );
% clim([-val,val]);
title('TB 2nd');
colorbar

diff_slice2=data(:,:,ix)-data2(:,:,ix);

val=0.5
figure
imagesc(diff_slice2);title('SB SB diff');
% clim([-val,val])
colorbar

min(diff_slice2,[],'all')
max(diff_slice2,[],'all')
min(data2,[],'all')
max(data2,[],'all')

ss=1