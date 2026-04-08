clear all
close all

root='/Volumes/ssd1/snapshots/Test_2500/';
% root='/Volumes/ssd1/snapshots/Test_2000/';
root='/Volumes/ssd1/snapshots/TEST_2514_no_check/'
fname=['TB_1st.raw'];
fname2=['TB_2nd.raw'];

root='/Volumes/ssd1/snapshots/TEST_source/'
fname=['SB_1st.raw'];
fname2=['SB_2nd.raw'];

root='/Volumes/ssd1/snapshots/TEST_accuracy/'
fname= ['TB_1st_1001.raw'];
fname2=['SB_1st_1000.raw'];

fname=fullfile(root,fname);
fname2=fullfile(root,fname2);

% add 8 points to each dimension
dims=[520,520,520];
ccnt=dims(1)*dims(2)*dims(3);
% 
% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data = reshape(A,dims);

fileID = fopen(fname2,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data2 = reshape(A,dims);
% 
data=permute(data, [3 2 1]);
data2=permute(data2, [3 2 1]);
close all
% 
ix=260;
% plotting
figure
imagesc( squeeze(data(:,:,ix)).'  );
% val=1.5;clim([-val,val]);
title(fname);
colorbar

figure
imagesc( squeeze(data2(:,:,ix)).' );
% clim([-val,val]);
title(fname2);
colorbar
% 
val=0.5;
diff_slice2=data(:,:,ix)-data2(:,:,ix);
figure
imagesc(diff_slice2);title('SB TB 1st diff');
colorbar

figure
imagesc(diff_slice2);title('SB TB 1st diff');
clim([0.4,val])
colorbar

% 
ss=1
% 
% diff_cube=data-data3;
% min(diff_cube,[],'all')
% max(diff_cube,[],'all')
% 
% figure
% imagesc(squeeze(diff_cube(:,ix,:)));title('SB TB 1st dif');
% clim([0.4,val])
% colorbar
% 
% figure
% imagesc(diff_cube(ix,:,:));title('SB TB 1st ');
% colorbar
% max(diff_slice4,[],'all')
% min(data2,[],'all')
% max(data2,[],'all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
