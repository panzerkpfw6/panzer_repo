clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% root='../../tb-sismos/SB_abc/data';
% root2='../../tb-sismos/TB_abc/data';
% fname= ['sismos.raw'];
% fname2=['sismos.raw'];

root='../../tb_modelling/SB_abc/data';
root2='../../tb_modelling/TB_abc/data';
fname= ['sismos.raw'];
fname2=['sismos.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname=fullfile(root,fname);
fname2=fullfile(root2,fname2);
%%%%%%%%%%%% xyz ordering
% dims_SB=[256,256,504];
% dims_TB=[256,256,505];
% dims_SB=[256,256,56];
% dims_TB=[256,256,57];
% dims_SB=[256,256,208];
% dims_TB=[256,256,209];
dims_SB=[256,512,504];
dims_TB=[256,512,505];
%%%%%%%%%%%% zxy ordering
% dims_SB=[504,256,512];
% dims_TB=[505,256,512];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ccnt_SB=dims_SB(1)*dims_SB(2)*dims_SB(3);
ccnt_TB=dims_TB(1)*dims_TB(2)*dims_TB(3);
s = dir(fname);      
filesize = s.bytes;
s2 = dir(fname2);         
filesize2 = s2.bytes;

% disp(strcat('theoretical file size= ',int2str(domain_size)));
disp(strcat('true file size1= ',int2str(filesize)));
disp(strcat('true file size2= ',int2str(filesize2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt_SB,'float32');
fclose(fileID);
dataSB = reshape(A,dims_SB);

fileID = fopen(fname2,'r');
A=fread(fileID,ccnt_TB,'float32');
fclose(fileID);
dataTB = reshape(A,dims_TB);
dataTB=dataTB(:,:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normalization
dataSB_normalized=dataSB./max(abs(dataSB),[],'all');
dataTB_normalized=dataTB./max(abs(dataTB),[],'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M,I] = max(dataSB,[],1);
[mxv,idx] = max(dataSB(:));
[r,c,p] = ind2sub(dims_SB,idx);
% dataSB(r,c,p)

min(dataSB,[],'all')
max(dataSB,[],'all')

min(dataTB,[],'all')
max(dataTB,[],'all')
% dataTB=permute(dataTB, [2 1 3]);
% dims_TB=[dims_TB(2),dims_TB(1),dims_TB(3)];

[M,I] = max(abs(dataTB),[],1);
[mxv,idx] = max(abs(dataTB(:)));
[r,c,p] = ind2sub(dims_TB,idx);
% dataTB(r,c,p)

% data_diff=dataSB-dataTB;
% min(data_diff,[],'all')
% max(data_diff,[],'all')
% max(data_diff);
% [M,I] = max(data_diff,[],1);
% [mxv,idx] = max(data_diff(:));
% [r,c,p] = ind2sub(size(data_diff),idx);
%%
project_dir='./';
font_sz=18;
dx=20;
dy=20;
dz=20;
dt=0.001;
x=(0:(dims_SB(1)-1))*dx/1000;
y=(0:(dims_SB(2)-1))*dy/1000;
z=(0:(dims_SB(3)-1))*dz/1000;
t=(0:(dims_SB(3)-1))*dt;

iiy=floor(dims_SB(2)/2);
iix=129;
% iix=10;
% iix=490;
iiy=200;
%%
fig=figure('units','normalized','outerposition',[0 0 0.8 0.46])  %0.8
val=1e-3;
ax1=subplot(1,2,1);
imagesc(ax1,y,t,squeeze(dataSB_normalized(iix,:,:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure');
% caxis([-val,val])
title('YZ shot gather,SB');
xlabel('Y, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');

ax2=subplot(1,2,2);
imagesc(ax2,y,t,squeeze(dataTB_normalized(iix,:,:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure');
% caxis([-val,val])
title('YZ shot gather,TB')
xlabel('Y, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_dir,'sismos.png'); 
delete(name); 
print(name,'-dpng',['-r' num2str(400)]);
% saveas(name)
ss=1
%%
fig=figure('units','normalized','outerposition',[0 0 0.8 0.46])  %0.8
val=1e-3;
ax1=subplot(1,2,1);
imagesc(ax1,x,t,squeeze(dataSB_normalized(:,iiy,:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure');
% caxis([-val,val])
title('XZ shot gather,SB');
xlabel('X, km');
ylabel('Time, sec');
set(gca,'FontSize',font_sz,'fontWeight','bold');

ax2=subplot(1,2,2);
imagesc(ax2,x,t,squeeze(dataTB_normalized(:,iiy,:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure');
% caxis([-val,val])
title('XZ shot gather,TB')
xlabel('X, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_dir,'sismos.png'); 
delete(name); 
print(name,'-dpng',['-r' num2str(400)]);
% saveas(name)
ss=1
%%
iiz=100;
figure
imagesc(squeeze(dataTB_normalized(:,:,iiz)).');
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure')
% caxis([-val,val])
title('XY shot gather,TB')
xlabel('X, km')
ylabel('Y, km')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_dir,'test.png');
delete(name); 
print(name,'-dpng',['-r' num2str(400)]);
%
figure
imagesc(squeeze(dataSB_normalized(:,:,iiz)).');
% imagesc(ax2,x,t,squeeze(dataTB(:,floor(dims_SB(2)/2),:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure')
% caxis([-val,val])
title('XY shot gather,SB')
xlabel('X, km')
ylabel('Y, km')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_dir,'test.png'); 
delete(name);
print(name,'-dpng',['-r' num2str(400)]);
%
% figure
% imagesc(squeeze(dataSB_normalized(:,iiz,:)).');
% % imagesc(ax2,x,t,squeeze(dataTB(:,floor(dims_SB(2)/2),:)).' );
% colormap('gray');
% hcb=colorbar;
% % hcb.Ticks=[];
% title(hcb,'P, pressure')
% % caxis([-val,val])
% title('XY shot gather,SB')
% xlabel('X, km')
% ylabel('Y, km')
% set(gca,'FontSize',font_sz,'fontWeight','bold');
% name=fullfile(project_dir,'test.png'); 
% delete(name);
% print(name,'-dpng',['-r' num2str(400)]);
%%
figure
imagesc(squeeze(dataSB(:,:,iiz)).'-squeeze(dataTB(:,:,iiz)).');
% imagesc(ax2,x,t,squeeze(dataTB(:,floor(dims_SB(2)/2),:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure')
% caxis([-val,val])
title('diff XY shot gather,TB')
xlabel('X, km')
ylabel('Y, km')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_dir,'test.png'); 
delete(name);
print(name,'-dpng',['-r' num2str(400)]);
%%
data_diff=dataSB-dataTB;
min(data_diff,[],'all')
max(data_diff,[],'all')

[mxv,idx] = max(data_diff(:));
[r,c,p] = ind2sub(dims_SB,idx);
data_diff(r,c,p)
disp([r,c,p]);
bb=1;

