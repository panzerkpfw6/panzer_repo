clear all
close all

root='/Volumes/ssd1/snapshots/Test_2500/';
% root='/Volumes/ssd1/snapshots/Test_2000/';
root='/Volumes/ssd1/snapshots/TEST_2514_no_check/'
root='/Volumes/ssd1/snapshots/TEST_source/'

fname=['TB_1st.raw'];
fname2=['TB_2nd.raw'];

root='/Volumes/ssd1/snapshots/TEST_source2/';
fname=['SB_1st_500.raw'];
fname2=['SB_2nd_500.raw'];

root='/Volumes/ssd1/snapshots/TEST_source3/';
fname=['SB_1st.raw'];
fname2=['SB_2nd.raw'];
% 
% root='/Volumes/ssd1/snapshots/TEST_source4/';
% root='/Users/pavelplotnitskii/Downloads/';
root='/Volumes/ssd1/snapshots/TEST_source5/'; %
% root='/Volumes/ssd1/snapshots/TEST_source6/'; 
root='/Volumes/ssd1/snapshots/TEST_long_int/'; 
fname=['TB_1st.raw']; fname2=['SB_1st.raw'];

root='/Volumes/ssd1/snapshots/TEST_source_6_CFL/'; 
fname=['SB_1st.raw']; fname2=['SB_2nd.raw'];

root='/Volumes/ssd1/snapshots/TEST_source_7_CFL/'; 

root='/Volumes/ssd1/snapshots/TEST_source_8_CFL/'; % large domain
fname=['SB_2nd.raw']; fname2=['SB_2nd.raw'];

fname=fullfile(root,fname);
fname2=fullfile(root,fname2);

% add 8 points to each dimension
dims=[520,520,520];
dims=[308,308,308];
dims=[1508,1508,1508];
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
ix=dims(1)/2;
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
% clim([0.4,val])
colorbar

figure;
plot(data(ix:end,ix,ix));hold on
plot(data2(ix:end,ix,ix));
legend('SB_1st','SB_2nd');

x1=data( ix:end,ix,ix);
x2=data2(ix:end,ix,ix);

figure;
% x2_tmp=circshift(x2,-1);
% plot( x1-x2_tmp );hold on

x2_tmp=circshift(x2,0);
plot( x1-x2_tmp );

% x2_tmp=circshift(x2,1);
% plot( x1-x2_tmp );
legend('1','2','3','4','5');
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
