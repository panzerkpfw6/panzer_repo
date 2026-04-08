clear all
close all

% dims=[128,256,512];
dims=[256,128,128];
dims=[676,676,201];

root='/media/plotnips/sdd1/Dropbox/Apps/download_to_laptop/salt_sh3_sb_rtm/'
root='/media/plotnips/sdd1/Dropbox/Apps/download_to_laptop/salt_sh3_sb_rtm/att4'
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/benchmark_512_3d.raw'];
fname=['/home/plotnips/Dropbox/PhD_proposal/work_with_david/Exawave_3_handover/simwave_export_to_ecrc_servers_/data/augmented_vel.raw'];
fname=['../rtm_munich/data/augmented_vel.raw'];
% fname=['../../stencil-main/data/used_vel_model.raw'];
fname=['../../data/velocity.raw'];
% fname=fullfile(root,'velocity.raw');
picture_dir='/media/plotnips/sdd1/Dropbox/Apps/Overleaf/2025_stencil_rtm/figures/rtm_salt/'

% some file size calculations
s = dir(fname);         
filesize = s.bytes ;
domain_size=dims(1)*dims(2)*dims(3)*4;
format longG
disp(strcat('true size=',int2str(filesize)))
disp(strcat('calculated domain size=',int2str(domain_size)))
data=read_snap(fname,'stencil',dims);   %velocity
[M,I] = max(data,[],1);
[mxv,idx] = max(abs(data(:)));
[r,c,p] = ind2sub(size(data),idx);
data1_max_loc=[r,c,p,mxv];
%% simwave indexing
dims2=dims;
fname2=['../rtm_munich/data/img.raw'];
% fname2=['../rtm_munich/data/backup2/img_1st.raw'];
% fname2=['../rtm_munich/data/ilm_only.raw'];
% fname2=['../rtm_munich/data/img_only.raw'];
% fname2=['../rtm_munich/data/img_2nd.raw'];
% fname2=['../rtm_munich/data/ilm_16455.raw'];
fname2=['../../data/img_1632.raw'];
% fname2=['../rtm_munich/data/snap_16449.raw'];
fname2=['../../data/ilm_only.raw'];     % img_only.raw
fname2=['../../data/img_only.raw'];     % img_only.raw
% fname2=['../../data/img_20472.raw'];     % img_only.raw
% fname2=['../../data/ilm_20468.raw'];
fname2=fullfile(root,'img_only.raw');
fname2=fullfile(root,'img.raw');
% fname2=fullfile(root,'ilm_only.raw');

s2=dir(fname2);filesize2=s2.bytes;disp(strcat('true size2=',int2str(filesize2)))
ccnt=dims2(1)*dims2(2)*dims2(3)*4;

rtm_data=read_snap(fname2,'stencil',dims2);   %rtm
[M,I] = max(rtm_data,[],1);
[mxv,idx] = max(abs(rtm_data(:)));
[r,c,p] = ind2sub(size(rtm_data),idx);
data2_max_loc=[r,c,p,mxv];
%% check
s=dir(fname);         
filesize1=s.bytes
s=dir(fname2);         
filesize2=s.bytes
dd=1
%%
% z1=25;
% data=data(:,:,z1:end);
% rtm_data=rtm_data(:,:,z1:end);
%% delete over-amplified top of rtm image
% data=data(:,:,20:end);
% rtm_data=rtm_data(:,:,20:end);
% dims=size(data);
% dims2=size(rtm_data);
%% coordinates
dx=25;dy=25;dz=25;
x=(0:(dims2(1)-1))*dx/1000;
y=(0:(dims2(2)-1))*dy/1000;
z=(0:(dims2(3)-1))*dz;
[x_,y_,z_] = meshgrid(x,y,z);
%% %% plotting
val1=1e-8; val2=1e-3;
val1=-1e-4; val2=-val1;
iy1=30;iy2=230;
ix=303;
%
figure
% imagesc( squeeze(rtm_data(128,:,:)).');
% imagesc( squeeze(rtm_data(303,:,:)).');
imagesc( squeeze(rtm_data(303,:,:)).',[val1,val2]);
% imagesc( squeeze(rtm_data(:,128,:)).',[-val,val] );
% imagesc( squeeze(rtm_data(:,133,:)).');
colormap('gray');
% clim([val1,val2])
colorbar
%%
% Clip amplitudes above a threshold (e.g., 95th percentile)
% threshold = prctile(rtm_data(:), 94); % 
threshold = prctile(rtm_data(:), 96); % my final choice for paper
% threshold=300
clipped_rtm = min(rtm_data, threshold); % Cap high values

% data2=clipped_rtm;
data2=rtm_data;
% % Optional: Apply a small Gaussian filter to smooth transitions
% data2 = imgaussfilt3(clipped_rtm, 1);

ix=413;
iy=300;
val1=-9.5*1e-0; % img_only.raw scaling
val1=-5.5*1e+2; % img_only.raw scaling
val1=-threshold
val2=-val1;

figure
subplot(1,2,1);
imagesc(y,z,squeeze(data(ix,:,:)).');
ax1=gca();
title(strcat('Velocity (m/sec), x=',num2str(ix)), 'Interpreter', 'latex' );
xlabel('Y, m');ylabel('Z, m');
colormap(ax1,'parula');
c = colorbar;

subplot(1,2,2);
imagesc(y,z,squeeze(data2(ix,:,:)).' ,[val1,val2] );
ax2=gca();
title(strcat('RTM image, x=',num2str(ix)));
xlabel('Y, m');ylabel('Z, m');
colormap(ax2,'gray');
colorbar
%
% val1=-1.0*1e-3; 
% val2=-val1;
figure
subplot(1,2,1);
imagesc(x,z, squeeze(data(:,iy,:)).');
ax1=gca();
title(strcat('$\mathbf{V_p}$ (m/sec), y=',num2str(iy)), 'Interpreter', 'latex','FontWeight', 'bold' );
xlabel('X, m');ylabel('Z, m');
colormap(ax1,'parula');
c = colorbar;

subplot(1,2,2);
imagesc(x,z, squeeze(data2(:,iy,:)).' ,[val1,val2] );
ax2=gca();
title(strcat('RTM image, y=',num2str(iy)));
xlabel('X, m');ylabel('Z, m');
colormap(ax2,'gray');
colorbar
%%
xslice = [1.000,8.000];    % location of y-z planes
yslice = [6.000,8.000,12.000,15.000];              % location of x-z plane
% zslice = [4500];         % location of x-y planes
% yslice = [];
zslice = [];

f1=figure('units','pixels','position',[2003 303 1605 558]);
h = slice(x_,y_,z_, data2, xslice, yslice, zslice);
set(h, 'FaceAlpha', 0.9)  % Add transparency to slices
ax = gca;  % Get the current axes
set(ax, 'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off')
set(ax, 'ZDir', 'reverse')
xlabel('X, km', 'FontWeight', 'bold');
ylabel('Y, km', 'FontWeight', 'bold');
zlabel('Z, m', 'FontWeight', 'bold');
grid off  % Remove grid lines
zlim([300 4000])


% set(ax, 'FontSize', 14)
% set(ax, 'XTickLabel', get(ax, 'XTickLabel'), 'FontSize', 12)
% set(ax, 'YTickLabel', get(ax, 'YTickLabel'), 'FontSize', 12)
% set(ax, 'ZTickLabel', get(ax, 'ZTickLabel'), 'FontSize', 12)
title('3D RTM image')

title(strcat('$\textbf{3D RTM image}$'),'Interpreter','latex','FontWeight','bold','FontSize', 19);
xlabel('$\textbf{X, km}$','Interpreter','latex','FontWeight','bold');
ylabel('$\textbf{Y, km}$','Interpreter','latex','FontWeight','bold');
zlabel('$\textbf{Z, m}$','Interpreter','latex','FontWeight','bold');
set(ax, 'FontSize', 19)

caxis([val1 val2])
c = colorbar;
c.Ticks = [];
colormap('gray')  % Use a colormap to highlight data variations

xlim([0 16.875])
ylim([0 16.875])

set(h,'edgecolor','none')
view(70.930049261083752,29.243382352941175)
% % Save the figure with specified DPI using print
print(f1, fullfile(picture_dir, 'rtm_vol'), '-dpng', '-r400');
%%
f1=figure('units','pixels','position',[2003 303 1200 600]);
imagesc(x,z, squeeze(data(:,iy,:)).');
ax1=gca();
title(strcat('$\mathbf{V_p (m/sec)},$','\textbf{ y=}',num2str((iy-1)*dy/1000,'%.1f'),'\textbf{ km}' ),'Interpreter','latex','FontWeight','bold','FontSize', 19);
xlabel('$\textbf{X, km}$','Interpreter','latex','FontWeight','bold');
ylabel('$\textbf{Z, m}$','Interpreter','latex','FontWeight','bold');
colormap(ax1,'parula');
set(ax1,'FontSize',19);
c = colorbar;
c.Ticks = [];
% Minimize white space around the axes
% Option 1: Tighten the axes margins using 'TightInset'
set(ax1, 'LooseInset', get(ax1, 'TightInset') + [0.02 0.02 0.02 0.02]); % Add small padding
% Adjust figure for saving
set(f1, 'PaperPositionMode', 'auto'); % Ensure the saved figure matches the screen
set(f1, 'InvertHardcopy', 'off'); % Preserve background color when saving
% print(f1, fullfile(picture_dir, 'vel_slice_'), '-dpng', '-r400');
exportgraphics(f1, fullfile(picture_dir,'vel_slice_.png'), 'Resolution', 400, 'BackgroundColor', 'white');
%%
f1 = figure('units', 'pixels', 'position', [2003 303 1200 600]);
imagesc(x, z, squeeze(data2(:, iy, :)).', [val1, val2]);
ax1 = gca();
title(strcat('\textbf{RTM image, y=}', num2str((iy-1)*dy/1000, '%.1f'), '\textbf{ km}'), ...
      'Interpreter', 'latex', 'FontWeight', 'bold', 'FontSize', 19);
xlabel('$\textbf{X, km}$', 'Interpreter', 'latex', 'FontWeight', 'bold');
ylabel('$\textbf{Z, m}$', 'Interpreter', 'latex', 'FontWeight', 'bold');
colormap(ax1,'gray');
set(ax1, 'FontSize', 19);
c = colorbar;
c.Ticks=[];
% Minimize white space around the axes
% Option 1: Tighten the axes margins using 'TightInset'
set(ax1, 'LooseInset', get(ax1, 'TightInset') + [0.02 0.02 0.02 0.02]); % Add small padding
% Adjust figure for saving
set(f1, 'PaperPositionMode', 'auto'); % Ensure the saved figure matches the screen
set(f1, 'InvertHardcopy', 'off'); % Preserve background color when saving
% print(f1, fullfile(picture_dir, 'rtm_slice_'), '-dpng', '-r400');
exportgraphics(f1, fullfile(picture_dir, 'rtm_slice_.png'), 'Resolution', 400, 'BackgroundColor', 'white');
%%
ss=1

function data=read_snap(fname,ordering,dims)
    ccnt=dims(1)*dims(2)*dims(3);
    fileID = fopen(fname,'r');
    A=fread(fileID,ccnt,'float32');
    fclose(fileID);

    s = dir(fname);
    filesize = s.bytes;

    if ordering=='stencil'
        nz=dims(3);
        data=nan(dims);
        tmp=1;
        for i=1:dims(1)    %x
            for j=1:dims(2)    %y
                data(i,j,:)=A(tmp:(tmp+nz-1));    %z
                tmp=tmp+nz;
            end
        end
    elseif ordering=='simwave'
        nx=dims(1);
        data=nan(dims);
        tmp=1;
        for k=1:dims(3)    %z
            for j=1:dims(2)    %y
                data(:,j,k)=A(tmp:(tmp+nx-1));    %x
                tmp=tmp+nx;
            end
        end
    end
end
