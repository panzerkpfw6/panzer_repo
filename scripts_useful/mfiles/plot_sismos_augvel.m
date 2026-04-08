clear all
close all
fname=['../simwave_export_to_ecrc_servers/data/sismos_52.raw'];
fname=['../simwave_export_to_ecrc_servers/data/sismos_53.raw'];
fname=['../simwave_export_to_ecrc_servers/data/sismos_54.raw'];
fname=['../simwave_export_to_ecrc_servers/data/sismos_3333.raw'];
dims=[51,51,100];
fname=['/Volumes/ssd1/SIMWAVE/simwave/seismic/sismos_6865.txt'];
dims=[11730,100];
dims=[11730,112];

fname=['/Volumes/ssd1/snapshots/sismos_1301.raw'];
dims=[11730,100];
dims=[11730,112];
dims=[2601,1000]


fname=['../../stencil-sismos/SB_abc/data/sismos_1302.raw'];
dims=[512,512,500]
dims=[676,676,1200]
dims=[676,676,4500]

fname=['../../stencil-sismos/SB_abc/data/sismos_1304.raw'];
dims=[200,200,2500]
dims=[512,512,200]
dims=[676,676,2800]
fname=['../../stencil-sismos/SB_abc/data/sismos.raw'];
dims=[676,676,3000]

% fname=['../../stencil-sismos/SB_abc/data/sismos_1301_simwave.raw'];
% dims=[512,512,1000]
% 
% fname=['../simwave_export_to_ecrc_servers_/data/sismos_131073.raw'];
% dims=[512,512,1500]

fname=['../../stencil-sismos/SB_abc/data/sismos_paper.raw'];
dims=[676,676,3000]


ccnt=dims(1)*dims(2)*dims(3);

% some file size calculations
s = dir(fname);         
filesize = s.bytes 
domain_size=ccnt*4;
format longG
disp(strcat('theoretical file size= ',int2str(domain_size)))
disp(strcat('true file size= ',int2str(filesize)))


% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
% A = fscanf(fileID,'%f',ccnt);
fclose(fileID);
data = reshape(A,dims);
data_normalized=data./max(abs(data),[],'all');
%%
% data2=permute(data, [3 1 2]);
% data2=zeros(dims);
% nlayer=dims(1)*dims(2);
% nt=dims(3);
% % A2=reshape(A,[nlayer,nt]);
% for k=1:dims(3)
%     disp(['k=',num2str(k)]);
%     for j=1:dims(2)
%         for i=1:dims(1)
% %             disp(['k=',num2str(k),',raw index=',num2str( (k-1)*nlayer+(j-1)*dims(2)+i )]);
%             data2(i,j,k)=A((k-1)*nlayer+(j-1)*dims(2)+i);
%         end
%     end
% end
% dims(1)*dims(2)
% %%
% figure
% j=100
% imagesc( squeeze(data2(:,j,:)).' );
% val=1e-12;
% title('2nd slice')
% colorbar
% xlabel('x')
% ylabel('t')
% 
% figure
% j=100
% imagesc( squeeze(data2(j,:,:)).' );
% val=1e-12;
% % clim([-val,val])
% title('1st slice')
% colorbar
% xlabel('x')
% ylabel('t')
% 
% figure
% j=100
% imagesc( squeeze(data2(:,:,100)).' );
% val=1e-12;
% % clim([-val,val])
% title('3rd slice')
% colorbar
% xlabel('x')
% ylabel('y')
%%
project_name='/home/plotnips/Dropbox/Apps/Overleaf/PASC_2024_stencil/figures';
font_sz=18;
dx=20;
dy=20;
dz=20;
dt=0.001;
x=(0:(dims(1)-1))*dx/1000;
y=(0:(dims(2)-1))*dy/1000;
z=(0:(dims(3)-1))*dz/1000;
t=(0:(dims(3)-1))*dt;

fig=figure('units','normalized','outerposition',[0 0 0.8 0.46])  %0.8
val=1e-3;
ax1=subplot(1,2,1);
imagesc(ax1,x,t,squeeze(data_normalized(:,floor(dims(2)/2),:)).' );
colormap('gray');
hcb=colorbar;hcb.Ticks=[];
title(hcb,'P, pressure');
caxis([-val,val])
title('XZ shot gather')
xlabel('Y, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');

ax2=subplot(1,2,2);
imagesc(ax2,y,t,squeeze(data_normalized(floor(dims(1)/2),:,:)).' );
colormap('gray');
hcb=colorbar;hcb.Ticks=[];
title(hcb,'P, pressure')
caxis([-val,val])
title('YZ shot gather')
xlabel('X, km')
ylabel('Time, sec')
set(gca,'FontSize',font_sz,'fontWeight','bold');
name=fullfile(project_name,'sismos.png'); 
delete(name); 
% print(name,'-dpng',['-r' num2str(400)]);
saveas(name)


figure
j=100
imagesc( squeeze(data_normalized(:,:,dims(3))).' );
colormap('gray');
colorbar
caxis([-val,val])
title('XY shot slice')
xlabel('x')
ylabel('y')
ss=1
%%
% for j=1:51
% for j=250:270
for j=90:110
    % plotting
    figure
    imagesc( squeeze(data(:,j,:)).' );
%     imagesc( squeeze(data(j,:,:)).' );
    val=1e-12;
    % clim([-val,val])
    title(j)
    colorbar
    xlabel('Y axis, grid cells')
    ylabel('T axis')
    w=waitforbuttonpress;
    close all
end

% figure
% imagesc( squeeze(data(:,25,:)).' );
% val=1e-12;
% % clim([-val,val])
% colorbar

min(data2,[],'all')
max(data2,[],'all')
min(data_normalized,[],'all')
max(data_normalized,[],'all')

ss=1