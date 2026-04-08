clear all
close all
pwd
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/marmousi_2301x512x751_xyz.raw'];
fname=['/Volumes/ssd1/velocity_models/marmousi_2301x512x751_xyz.raw'];
dims=[2301,512,751];
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

slice1=squeeze(data(:,310,:)).';
slice2=squeeze(data(:,510,:)).';
diff=slice1-slice2;

figure
imagesc( slice1 );
colorbar

font_sz=18;
figure
imagesc( slice2 );
xlabel('X axis');
ylabel('Z axis');
colorbar
set(gca,'FontSize',font_sz,'fontWeight','bold');
name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig/xz_marm.png'; delete(fname)
print(name,'-dpng',['-r' num2str(400)]);

figure
imagesc( squeeze(data(1100,:,:)).' );
xlabel('Y axis');
ylabel('Z axis');
colorbar
set(gca,'FontSize',font_sz,'fontWeight','bold');
name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig/yz_marm_3d.png'; delete(fname)
print(name,'-dpng',['-r' num2str(400)]);

figure
imagesc(diff);
colorbar

ss=1;
