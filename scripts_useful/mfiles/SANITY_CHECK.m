clear all
close all

root='/Volumes/ssd1/snapshots/Test_2500/';
% root='/Volumes/ssd1/snapshots/Test_2000/';
root='/Volumes/ssd1/snapshots/TEST_2514_no_check/'

fname=['SB_1st.raw'];
fname2=['SB_1st_abc.raw'];
fname3=['TB_1st.raw'];
fname4=['TB_1st_abc.raw'];
fname5=['SB_2nd.raw'];

% fname=fullfile(root,fname);
% fname2=fullfile(root,fname2);
% fname3=fullfile(root,fname3);
% fname4=fullfile(root,fname4);
% fname5=fullfile(root,fname5);
% 
% % add 8 points to each dimension
% dims=[520,520,520];
% ccnt=dims(1)*dims(2)*dims(3);
% 
% % read file according to dimensions
% fileID = fopen(fname,'r');
% A=fread(fileID,ccnt,'float32');
% fclose(fileID);
% data = reshape(A,dims);
% 
% fileID = fopen(fname2,'r');
% A=fread(fileID,ccnt,'float32');
% fclose(fileID);
% data2 = reshape(A,dims);
% 
% fileID = fopen(fname3,'r');
% A=fread(fileID,ccnt,'float32');
% fclose(fileID);
% data3 = reshape(A,dims);
% 
% fileID = fopen(fname4,'r');
% A=fread(fileID,ccnt,'float32');
% fclose(fileID);
% data4 = reshape(A,dims);
% 
% data=permute(data, [3 2 1]);
% data2=permute(data2, [3 2 1]);
% data3=permute(data3, [3 2 1]);
% data4=permute(data4, [3 2 1]);
% close all
% 
% ix=260;
% % plotting
% figure
% imagesc( squeeze(data(:,:,ix)).'  );
% % val=1.5;clim([-val,val]);
% title(fname);
% colorbar
% 
% figure
% imagesc( squeeze(data2(:,:,ix)).' );
% % clim([-val,val]);
% title(fname2);
% colorbar
% 
% figure
% imagesc( squeeze(data3(:,:,ix)).' );
% % clim([-val,val]);
% title(fname3);
% colorbar
% 
% figure
% imagesc( squeeze(data4(:,:,ix)).' );
% % clim([-val,val]);
% title(fname4);
% colorbar
% 
% diff_slice2=data(:,:,ix)-data2(:,:,ix);
% diff_slice3=data(:,:,ix)-data3(:,:,ix);
% diff_slice4=data(:,:,ix)-data4(:,:,ix);
% diff_slice5=data2(:,:,ix)-data4(:,:,ix);
% 
% figure
% imagesc(diff_slice3);title('SB TB 1st diff');
% colorbar
% figure
% imagesc(diff_slice5);title('SB TB 1st abc diff');
% colorbar
% 
% val=0.5
% figure
% imagesc(diff_slice3);title('SB TB 1st diff,[0.4,0.5]');
% clim([0.4,val])
% colorbar
% figure
% imagesc(diff_slice5);title('SB TB 1st abc diff,[0.4,0.5]');
% clim([0.4,val])
% colorbar
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
% 
% 
% min(diff_slice2,[],'all')
% max(diff_slice2,[],'all')
% min(diff_slice3,[],'all')
% max(diff_slice3,[],'all')
% min(diff_slice4,[],'all')
% max(diff_slice4,[],'all')
% 
% min(data2,[],'all')
% max(data2,[],'all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='/Volumes/ssd1/snapshots/Test_2500/'
root='/Volumes/ssd1/snapshots/Test_2000/';

fname= ['SB_2nd.raw'];
fname2=['SB_2nd_abc.raw'];
fname3=['TB_2nd.raw'];
fname4=['TB_2nd_abc.raw'];
fname5=['SB_1st.raw'];

fname=fullfile(root,fname);
fname2=fullfile(root,fname2);
fname3=fullfile(root,fname3);
fname4=fullfile(root,fname4);
fname5=fullfile(root,fname5);

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

fileID = fopen(fname3,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data3 = reshape(A,dims);

fileID = fopen(fname4,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data4 = reshape(A,dims);

data=permute(data, [3 2 1]);
data2=permute(data2, [3 2 1]);
data3=permute(data3, [3 2 1]);
data4=permute(data4, [3 2 1]);
close all

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

figure
imagesc( squeeze(data3(:,:,ix)).' );
% clim([-val,val]);
title(fname3);
colorbar

figure
imagesc( squeeze(data4(:,:,ix)).' );
% clim([-val,val]);
title(fname4);
colorbar

diff_slice2=data(:,:,ix)-data2(:,:,ix);
diff_slice3=data(:,:,ix)-data3(:,:,ix);
diff_slice4=data(:,:,ix)-data4(:,:,ix);
diff_slice5=data2(:,:,ix)-data4(:,:,ix);

figure
imagesc(diff_slice3);title('SB TB 1st diff');
colorbar
figure
imagesc(diff_slice5);title('SB TB 1st abc diff');
colorbar

val=0.5
figure
imagesc(diff_slice3);title('SB TB 1st diff,[0.4,0.5]');
clim([0.4,val])
colorbar
figure
imagesc(diff_slice5);title('SB TB 1st abc diff,[0.4,0.5]');
clim([0.4,val])
colorbar

diff_cube=data-data3;
min(diff_cube,[],'all')
max(diff_cube,[],'all')

figure
imagesc(squeeze(diff_cube(:,ix,:)));title('SB TB 1st dif');
clim([0.4,val])
colorbar

figure
imagesc(diff_cube(ix,:,:));title('SB TB 1st ');
colorbar


min(diff_slice2,[],'all')
max(diff_slice2,[],'all')
min(diff_slice3,[],'all')
max(diff_slice3,[],'all')
min(diff_slice4,[],'all')
max(diff_slice4,[],'all')

min(data2,[],'all')
max(data2,[],'all')
ss=1
