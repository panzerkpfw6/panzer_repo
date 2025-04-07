clear all
close all

% dims=[128,256,512];
dims=[256,128,128];
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/benchmark_512_3d.raw'];
fname=['/home/plotnips/Dropbox/PhD_proposal/work_with_david/Exawave_3_handover/simwave_export_to_ecrc_servers_/data/augmented_vel.raw'];
fname=['../rtm_munich/data/augmented_vel.raw'];
% fname=['../../stencil-main/data/used_vel_model.raw'];
fname=['../../data/velocity.raw'];

% some file size calculations
s = dir(fname);         
filesize = s.bytes ;
domain_size=2301*512*751*4;
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
% fname2=['../../data/img_only.raw'];     % img_only.raw
fname2=['../../data/img_1636.raw'];

s2=dir(fname2);filesize2=s2.bytes;disp(strcat('true size2=',int2str(filesize2)))
ccnt=dims2(1)*dims2(2)*dims2(3);

data2=read_snap(fname2,'stencil',dims2);   %rtm
[M,I] = max(data2,[],1);
[mxv,idx] = max(abs(data2(:)));
[r,c,p] = ind2sub(size(data2),idx);
data2_max_loc=[r,c,p,mxv];
%% %% plotting
val1=1e-8; val2=1e-3;
val1=-1e-1; val2=-val1;
iy1=30;iy2=230;
%
figure
imagesc( squeeze(data2(128,:,:)).');
% imagesc( squeeze(data2(:,128,:)).',[-val,val] );
% imagesc( squeeze(data2(:,133,:)).');
colormap('gray');
clim([val1,val2])
colorbar
%%
figure
% imagesc( squeeze(data2(6,:,:)).',[-val,val]);
imagesc( squeeze(data2(6,:,:)).');
colormap('gray');
colorbar
ss=1
figure
imagesc( squeeze(data(64,iy1:iy2,:)).',[-val,val]);
colormap('gray');
colorbar

figure
subplot(1,2,1);ax1=gca();
imagesc( squeeze(data(64,iy1:iy2,:)).');
colorbar

subplot(1,2,2);ax2=gca();
imagesc( squeeze(data2(64,iy1:iy2,:)).',[-val,val]);
colorbar
title(ax1,'Velocity (m/sec)');
title(ax2,'RTM image');
colormap(ax1,'parula');
colormap(ax2,'gray');
xlabel(ax1,'Y');ylabel(ax1,'Y');
xlabel(ax2,'Y');ylabel(ax2,'Y');


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
