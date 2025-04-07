clear all
close all

dims=[128,256,512];
dims=[128,256,128];
% dims=[128+8,256+8,512+8];
fname=['/Volumes/ssd1/SIMWAVE/simwave/models/benchmark_512_3d.raw'];
fname=['/home/plotnips/Dropbox/PhD_proposal/work_with_david/Exawave_3_handover/simwave_export_to_ecrc_servers_/data/augmented_vel.raw'];
fname=['../rtm_munich/data/augmented_vel.raw'];
% fname=['../../stencil-main/data/used_vel_model.raw'];

dims=[676,676,201];
fname=['../stencil-rtm/data/augmented_vel.raw'];
%%
dims=[256,256,256];
fname=['../../data/velocity.raw'];
fname=['../../data/density.raw'];
fname=['../../data/coef_1st.raw'];
% fname=['../../data/coef_2nd.raw'];
%%
dims=[676,676,201];
fname=['../../data/velocity.raw'];
%%
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
data_max_loc=[r,c,p,mxv];
%% %% plotting
val=1e+3;
iy1=30;iy2=230;
%
figure
ax1=gca();
% imagesc( squeeze(data(:,310,:)).');
imagesc( squeeze(data(310,:,:)).');
colorbar
title(ax1,'Velocity (m/sec)');
colormap(ax1,'parula');
xlabel(ax1,'Y');ylabel(ax1,'Y');

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
