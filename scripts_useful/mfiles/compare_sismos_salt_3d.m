clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='../../tb-sismos/SB_abc/data';
root2='../../tb-sismos/TB_abc/data';
fname= ['sismos.raw'];
fname2=['sismos.raw'];
fname3=['used_vel_model.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname=fullfile(root,fname);
fname2=fullfile(root2,fname2);
fname_vel=fullfile(root,fname3);

dims_vel=[640,640,201];
dims_SB=[256,256,504];
dims_TB=[256,256,505];
% dims_SB=[256,256,56];
% dims_TB=[256,256,57];
% dims_SB=[256,256,208];
% dims_TB=[256,256,209];
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read velocity model
fileID = fopen(fname_vel,'r');
A=fread(fileID,ccnt_SB,'float32');
fclose(fileID);
vel= reshape(A,dims_vel);

figure
imagesc( squeeze(vel(:,:,300)).'  );
% caxis([-val,val]);
title(fname3);
colorbar
figure
imagesc( squeeze(vel(300,:,:)).'  );
% caxis([-val,val]);
title(fname3);
colorbar
figure
imagesc( squeeze(vel(:,300,:)).'  );
% caxis([-val,val]);
title(fname3);
colorbar
%%
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
% dataTB=permute(dataTB, [2 1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normalization
dataSB_normalized=dataSB./max(abs(dataSB),[],'all');
dataTB_normalized=dataTB./max(abs(dataTB),[],'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M,I] = max(dataSB,[],1);
[mxv,idx] = max(dataSB(:));
[r,c,p] = ind2sub(size(dataSB),idx);
dataSB(r,c,p)

min(dataSB,[],'all')
max(dataSB,[],'all')

min(dataTB,[],'all')
max(dataTB,[],'all')

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

fig=figure('units','normalized','outerposition',[0 0 0.8 0.46])  %0.8
val=1e-3;
ax1=subplot(1,2,1);
imagesc(ax1,x,t,squeeze(dataSB_normalized(:,floor(dims_SB(2)/2),:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure');
% caxis([-val,val])
title('XZ shot gather,SB')
xlabel('Y, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');

ax2=subplot(1,2,2);
imagesc(ax2,x,t,squeeze(dataTB_normalized(:,floor(dims_SB(2)/2),:)).' );
% imagesc(ax2,x,t,squeeze(dataTB(:,floor(dims_SB(2)/2),:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure')
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
iiz=200;
%%
figure
imagesc(squeeze(dataTB_normalized(:,:,iiz)).');
% imagesc(ax2,x,t,squeeze(dataTB(:,floor(dims_SB(2)/2),:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure')
% caxis([-val,val])
title('XZ shot gather,TB')
xlabel('X, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_dir,'test.png');
delete(name); 
print(name,'-dpng',['-r' num2str(400)]);
%%
figure
imagesc(squeeze(dataSB_normalized(:,:,iiz)).');
% imagesc(ax2,x,t,squeeze(dataTB(:,floor(dims_SB(2)/2),:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure')
% caxis([-val,val])
title('XZ shot gather,TB')
xlabel('X, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_dir,'test.png'); 
delete(name);
print(name,'-dpng',['-r' num2str(400)]);
%%
figure
imagesc(squeeze(dataSB_normalized(:,:,iiz)).'-squeeze(dataTB_normalized(:,:,iiz)).');
% imagesc(ax2,x,t,squeeze(dataTB(:,floor(dims_SB(2)/2),:)).' );
colormap('gray');
hcb=colorbar;
% hcb.Ticks=[];
title(hcb,'P, pressure')
% caxis([-val,val])
title('diff XY shot gather,TB')
xlabel('X, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_dir,'test.png'); 
delete(name);
print(name,'-dpng',['-r' num2str(400)]);

bb=1;

