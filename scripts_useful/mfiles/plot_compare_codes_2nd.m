clear all
close all
19140-4*136
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1500*1500*0.001*0.001/10/10  % CFL stencil. should be <1
1500*1500*0.001*0.001/0.0025/0.0025 % simwave stencil. should be <1
1500*0.001/10  % CFL stencil. should be <1
c=1500;hh=10;dt=0.001;
3*c^2*dt^2/(hh^2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 512 512 512 size
% read SB stencil snapshot
root='../../stencil-main';
fname= ['SB_1st_abc.raw'];
% fname= ['SB_2nd_abc.raw'];
% read SB simwave snapshot
root2='../simwave_export_to_ecrc_servers_';
root2='../rtm_munich';
fname2=['snapshot_504'];
fname2=['TB_lastshot_u0'];
root2='../rtm_munich';
fname2=['snapshot_504'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname=fullfile(root,fname);
fname2=fullfile(root2,fname2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add 8 points to each dimension
% dims1=[264,264,264];
% dims1=[520,520,520];
dims1=[128+8,256+8,512+8];
dims2=[128+8,256+8,512+8];
% dims2=[512+8,256+8,128+8];
% dims2=dims1;
ccnt=dims1(1)*dims1(2)*dims1(3);
%%%%%%%%%%%%%%%%%%%%%%%%% check sizes
s = dir(fname);
filesize = s.bytes;
s2 = dir(fname2);
filesize2 = s2.bytes;
disp(filesize)
disp(filesize2)
%%%%%%%%%%%%%%%%%%%%%%%%% read file according to dimensions and reshape
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);

fileID = fopen(fname2,'r');
A2=fread(fileID,ccnt,'float32');
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%% Alternative reshaping, data
nz=dims1(3);
data=nan(dims1);
tmp=1;
for i=1:dims1(1)    %x
    for j=1:dims1(2)    %y
        data(i,j,:)=A(tmp:(tmp+nz-1));    %z
        tmp=tmp+nz;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%% Alternative reshaping, data2
nx=dims2(1);
data2=nan(dims2);
tmp=1;
for k=1:dims2(3)    %z
    for j=1:dims2(2)    %y
        data2(:,j,k)=A2(tmp:(tmp+nx-1));    %x
        tmp=tmp+nx;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%% reshape data
% data = reshape(A,dims1);
% data2= reshape(A2,dims2);
%%%%%%%%%%%%%%%%%%%%%%%%% permute data
% data=permute(data, [3 2 1]);
% not permute for now
% data2=permute(data2, [3 2 1]);
% data2=permute(data2,[2 1 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%
data_diff=data-data2;
1+1e-11-1
C =1+ 0.1.^(+12)-1


min(data_diff,[],'all')
max(data_diff,[],'all')
abs_diff=max(abs(data_diff),[],'all')
RMS_val = rms(data_diff,"all")

close all

[M,I] = max(data,[],1);
[mxv,idx] = max(abs(data(:)));
[r,c,p] = ind2sub(size(data),idx);
data(r,c,p)
data1_max_loc=[r,c,p,mxv];

[M,I] = max(data2,[],1);
[mxv,idx] = max(abs(data2(:)));
[r2,c2,p2] = ind2sub(size(data2),idx);
data2_max_loc=[r2,c2,p2,mxv];

% max(abs(data_diff(:)))
% [M,I] = max(data_diff,[],1);
% [mxv,idx] = max(data_diff(:));
% [r,c,p] = ind2sub(size(data_diff),idx);
% data_diff(r,c,p)

% ix=dims(3)/2;
ix=256;
val=1.5e-4;

% figure
% imagesc( squeeze(data_diff(:,:,ix)).' );
% title('diff 2nd SB TB');
% colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1e-3;
a=10;
r2=68;c2=132;p2=260;
figure
subplot(1,3,1)
imagesc( squeeze(data2(r2,:,:)).' );    %p2
% imagesc( squeeze(data2(r2,:,:)).' ,[-a,a] );    %p2
title(strcat('rtm,x=',num2str(r2)) );
xlabel('Y');ylabel('Z');
colorbar

subplot(1,3,2)
imagesc( squeeze(data2(:,c2,:)).' );    %p2
% imagesc( squeeze(data2(:,c2,:)).' ,[-a,a] );    %p2
title(strcat('rtm,y=',num2str(c2)));
xlabel('X');ylabel('Z');
colorbar

subplot(1,3,3)
imagesc( squeeze(data2(:,:,p2)).');    %p2
% imagesc( squeeze(data2(:,:,p2)).' ,[-a,a] );    %p2
title(strcat('rtm,z=',num2str(p2)));
xlabel('X');ylabel('Y');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=68;c=132;p=260;
figure
subplot(1,3,1)
imagesc( squeeze(data(r+1,:,:)).'  );
title(strcat('stencil,x=',num2str(r)));
xlabel('Y');ylabel('Z');
colorbar

subplot(1,3,2)
imagesc( squeeze(data(:,c+1,:)).'  );    %p
title(strcat('stencil,y=',num2str(c)));
xlabel('X');ylabel('Z');
colorbar

subplot(1,3,3)
imagesc( squeeze(data(:,:,p+1)).'  );    %p
title(strcat('stencil,z=',num2str(p)));
xlabel('X');ylabel('Y');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
title('Difference,\n second order solution')
subplot(1,3,1)
imagesc( squeeze(data_diff(r2,:,:)).' );    %p2
title(strcat('diff,x=',num2str(r2)) );
xlabel('Y');ylabel('Z');
colorbar

subplot(1,3,2)
imagesc( squeeze(data_diff(:,c2,:)).' );    %p2
title(strcat('diff,y=',num2str(c2)));
xlabel('X');ylabel('Z');
colorbar

subplot(1,3,3)
imagesc( squeeze(data_diff(:,:,p2)).' );    %p2
title(strcat('diff,z=',num2str(p2)));
xlabel('X');ylabel('Y');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_diff=max(abs(data_diff),[],'all')
RMS_val = rms(data_diff,"all")
ss=1