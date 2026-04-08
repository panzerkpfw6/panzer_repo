clear all
close all

root='/Volumes/ssd1/snapshots';
folder_name=['SB_1st_marm_3000_it'];
folder_name=['SB_1st_marm_abc_3000_it'];
% fname=['SB_1st_marm_abc_6000_it'];
dims=[2301,512,751];
%%%%%%%%%%%%%%%%%
folder_name=['SB_1st_salt_1000_nt'];
folder_name=['SB_1st_salt_nt_1000'];
folder_name=['SB_1st_salt_nt_2000'];
folder_name=['SB_1st_salt_nt_1500'];
folder_name=['SB_1st_salt3d_nt_1000'];
folder_name=['SB_1st_salt_nt1000_isz2'];
folder_name=['SB_1st_salt_nt1500_isz5'];
folder_name=['SB_1st'];
folder_name=['SB_1st_salt_nt1200_isz2'];
fname=strcat(folder_name,'.raw')
dims=[676,676,201];
fname=fullfile(root,fname);

% fname=['../../stencil-sismos/SB_abc/SB_1st_abc.raw']
% dims=[512,512,512];
%%%%%%%%%%%%%%%%%
% add 8 points to each dimension
dims=dims+8;
ccnt=dims(1)*dims(2)*dims(3);
%%
dx=20;
dy=20;
dz=20;
x=(0:(dims(1)-1))*dx/1000;
y=(0:(dims(2)-1))*dy/1000;
z=(0:(dims(3)-1))*dz/1000;
% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_orig= reshape(A,[dims(3),dims(2),dims(1)]);
data=permute(data_orig, [3 2 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig';
% folder_name='marm_3000_nt';
% folder_name='marm_abc_3000_nt';
% % folder_name='marm_abc_6000_nt';
% folder_name='salt_abc_2000_nt';
% folder_name='salt_abc_1500_nt';
% folder_name='salt_abc_1000_nt';
mkdir(fullfile(project_name,folder_name));

font_sz=18;
val=0.5;
close all
%%
% snap=squeeze(data(200,:,:)).';
% snap=snap./max(abs(snap),[],'all');
% figure; 
% imagesc(y,z,snap );
% xlabel('Y, km');
% ylabel('Z, km');
% % clim([-val,val]);
% colormap('gray');
% set(gca,'FontSize',font_sz,'fontWeight','bold')
% colorbar 
% pic_name='yz_snap.png';
% name=fullfile(project_name,folder_name,pic_name);
% delete(name)
% print(name,'-dpng',['-r' num2str(400)]);
% 
% snap=squeeze(data(:,100,:)).';
% snap=snap./max(abs(snap),[],'all');
% figure; imagesc( x,z,snap );
% % clim([-val,val]);
% colormap('gray');
% xlabel('X, km');
% ylabel('Z, km');
% set(gca,'FontSize',font_sz,'fontWeight','bold')
% colorbar
% pic_name='xz_snap.png';
% name=fullfile(project_name,folder_name,pic_name);
% delete(name)
% % print(name,'-dpng',['-r' num2str(400)]);
% 
% snap=squeeze(data(:,:,5)).';
% snap=snap./max(abs(snap),[],'all');
% figure; imagesc( x,y,snap );
% colormap('gray');
% xlabel('X, km');
% ylabel('Y, km');
% set(gca,'FontSize',font_sz,'fontWeight','bold')
% colorbar
% pic_name='xy_snap.png';
% name=fullfile(project_name,folder_name,pic_name);
% delete(name)
% % print(name,'-dpng',['-r' num2str(400)]);
% 
% snap=squeeze(data(:,:,100)).';
% snap=snap./max(abs(snap),[],'all');
% figure; imagesc( x,y,snap );
% colormap('gray');
% xlabel('X, km');
% ylabel('Y, km');
% set(gca,'FontSize',font_sz,'fontWeight','bold')
% colorbar
% pic_name='xy_snap_100.png';
% name=fullfile(project_name,folder_name,pic_name);
% delete(name)
% % print(name,'-dpng',['-r' num2str(400)]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
snap=squeeze(data(200,:,:)).';
snap=snap./max(abs(snap),[],'all');
slice_yz=snap;
snap=squeeze(data(:,100,:)).';
snap=snap./max(abs(snap),[],'all');
slice_xz=snap;
snap=squeeze(data(:,:,100)).';
snap=snap./max(abs(snap),[],'all');
slice_xy_100=snap;

figure('units','normalized','outerposition',[0 0 1.0 0.36])
ax1=subplot(1,3,1);
imagesc(ax1,x,z, slice_xz );
colormap('gray');
xlabel('X, km');
ylabel('Z, km');
title('XZ slice at Y=2 km');
hcb=colorbar
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString,'FontSize',font_sz);
set(gca,'FontSize',font_sz,'fontWeight','bold');

ax2=subplot(1,3,2);
imagesc(ax2,y,z, slice_yz );
colormap('gray');
xlabel('Y, km');
ylabel('Z, km');
title('YZ slice at X=4 km');
hcb=colorbar
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString,'FontSize',font_sz);
set(gca,'FontSize',font_sz,'fontWeight','bold');

ax3=subplot(1,3,3);
imagesc(ax3,x,y,slice_xy_100 );
colormap('gray');
xlabel('X, km');
ylabel('Y, km');
title('XY slice at Z=2 km')
hcb=colorbar
colorTitleHandle = get(hcb,'Title');
titleString = 'V, km/sec';
set(colorTitleHandle ,'String',titleString,'FontSize',font_sz);
set(gca,'FontSize',font_sz,'fontWeight','bold');
name='/Users/pavelplotnitskii/Library/CloudStorage/Dropbox/Apps/Overleaf/seg_23_alkhobar_workshop/Fig/snapshot2.png'; delete(name); print(name,'-dpng',['-r' num2str(400)]);
ss=1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some file size calculations
s = dir(fname);         
filesize = s.bytes 
domain_size=ccnt*4;
format longG
disp(domain_size)
disp(filesize)
pp=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% code backup%%%%%%%%%%%%%%%%%%%%%%%%
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
% figure; imagesc( squeeze(data(:,354,:)) );colorbar
% figure; imagesc( squeeze(data(:,90,:)) );colorbar
% figure; imagesc( squeeze(data(:,:,400)) );colorbar
% figure; imagesc( squeeze(data(254,:,:)) );colorbar
% figure; imagesc( squeeze(data(:,400,:)) );colorbar
% figure; imagesc( squeeze(data(:,330,:)) );colorbar
% figure; imagesc( squeeze(data(300,:,:)) );colorbar


% data=permute(data_orig, [1 3 2]);
% data=permute(data_orig, [2 3 1]);
% data=permute(data_orig, [2 1 3]);
% data=permute(data_orig, [3 2 1]);
% data=permute(data_orig, [3 1 2]);

% data= reshape(A,[dims(1),dims(2),dims(3)]);
% data= reshape(A,[dims(1),dims(3),dims(2)]);
% data= reshape(A,[dims(2),dims(3),dims(1)]);
% data= reshape(A,[dims(2),dims(1),dims(3)]);
% data= reshape(A,[dims(3),dims(2),dims(1)]);
% data= reshape(A,[dims(3),dims(1),dims(2)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% searching for data min and max 
% [M,I]=max(data,[],'all');
% [i1,i2,i3] = ind2sub(size(data),I);
% 
% figure;imagesc( squeeze(data(:,i2,:)).' );colorbar
% 
% figure;imagesc( squeeze(data(i1,:,:)).' );colorbar
% 
% figure;imagesc( squeeze(data(:,:,i3)).' );colorbar
% dd=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% code backup%%%%%%%%%%%%%%%%%%%%%%%%