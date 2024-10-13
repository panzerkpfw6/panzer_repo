clear all
close all
pwd
fname=['/Volumes/ssd1/SIMWAVE/simwave_mlbs/simwave/models/marmousi_2301x512x751_xyz.raw'];
fname=['/Volumes/ssd1/velocity_models/salt3d_676x676x201_zyx.raw'];
fname=['../../velocity_models/salt3d_676x676x201_xyz.raw'];
dims=[676,676,201];
ccnt=dims(1)*dims(2)*dims(3);

% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
A=A/1000;
%%
data= reshape(A,dims);
%%
dx=20;
dy=20;
dz=20;
x=(0:(dims(1)-1))*dx/1000;
y=(0:(dims(2)-1))*dy/1000;
z=(0:(dims(3)-1))*dz/1000;
%% record reshaped model
% data= reshape(A,[dims(3),dims(2),dims(1)]);
% data=permute(data, [3 2 1]);
% fname2=['/Volumes/ssd1/velocity_models/salt3d_676x676x201_xyz.raw'];
% fileID = fopen(fname2,'w');
% fwrite(fileID,data,'float32');
% fclose(fileID);
% ss=1;
%%
isx=200;
isy=100;
isz=2;
floor(dims(1)/2)
slice_yz=squeeze(data(floor(dims(1)/2),:,:)).';
slice_xz=squeeze(data(:,floor(dims(2)/2),:)).';
slice_xy=squeeze(data(:,:,10)).';
slice_xy_100=squeeze(data(:,:,100)).';


font_sz=18;
figure
imagesc(x,z, slice_xz );
xlabel('X, km');
ylabel('Z, km');
hcb=colorbar
set(gca,'FontSize',font_sz,'fontWeight','bold');
% name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig/xz_salt.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);
% name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/hpc_asia_2024/Fig/xz_salt.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);

figure
imagesc(y,z, slice_yz );
xlabel('Y, km');
ylabel('Z, km');
colorbar
set(gca,'FontSize',font_sz,'fontWeight','bold');
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString);
% name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig/yz_salt.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);
% name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/hpc_asia_2024/Fig/yz_salt.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);


figure
imagesc(x,y,slice_xy );
% imagesc( squeeze(data(:,:,100)).' );
xlabel('X, km');
ylabel('Y, km');
hcb=colorbar
set(gca,'FontSize',font_sz,'fontWeight','bold');
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString);
% name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig/xy_salt.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);

figure
imagesc(x,y,slice_xy_100 );
xlabel('X, km');
ylabel('Y, km');
hcb=colorbar
set(gca,'FontSize',font_sz,'fontWeight','bold');
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString);
% name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig/xy_salt_z100.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);
%%
font_sz=17;
figure('units','normalized','outerposition',[0 0 1.0 0.36])
ax1=subplot(1,3,1);
imagesc(ax1,x,z, slice_xz );
xlabel('X, km');
ylabel('Z, km');
title('XZ slice at Y=6.76 km.');
hcb=colorbar
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString,'FontSize',font_sz);
set(gca,'FontSize',font_sz,'fontWeight','bold');

ax2=subplot(1,3,2);
imagesc(ax2,y,z, slice_yz );
xlabel('Y, km');
ylabel('Z, km');
title('YZ slice at X=6.76 km.');
hcb=colorbar
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString,'FontSize',font_sz);
set(gca,'FontSize',font_sz,'fontWeight','bold');

ax3=subplot(1,3,3);
imagesc(ax3,x,y,slice_xy_100 );
xlabel('X, km');
ylabel('Y, km');
title('XY slice at Z=2 km.')
hcb=colorbar
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString,'FontSize',font_sz);
set(gca,'FontSize',font_sz,'fontWeight','bold');
% name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/PASC_2024_stencil/Fig/salt_vel3.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);
name='/home/plotnips/Dropbox/Apps/Overleaf/PASC_2024_stencil/Fig/salt_vel3.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);
%%
ss=1

% format longG
% s = dir(fname);
% filesize = s.bytes 
% disp(strcat('theoretical file size= ',int2str(ccnt*4)))
% disp(strcat('true file size= ',int2str(filesize)))
