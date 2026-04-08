clear all
close all
19140-4*136
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1500*1500*0.001*0.001/10/10  % CFL stencil. should be <1
1500*1500*0.001*0.001/0.0025/0.0025 % simwave stencil. should be <1
1500*0.001/10  % CFL stencil. should be <1
c=1500;hh=10;dt=0.001;
3*c^2*dt^2/(hh^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='../stencil-rtm-dev';
fname=['snapshot_400'];
fname=['snapshot_SB1st_400'];
fname=['snapshot_SB2nd_400'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname=fullfile(root,fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add 8 points to each dimension

dims1=[676+8,676+8,201+8];
dims1=[128+8,256+8,512+8];
ccnt=dims1(1)*dims1(2)*dims1(3);
%%%%%%%%%%%%%%%%%%%%%%%%% check sizes
s = dir(fname);
filesize = s.bytes;
disp(filesize)
%%%%%%%%%%%%%%%%%%%%%%%%% read file according to dimensions and reshape
data=read_snap(fname,'simwave',dims1);   %SB    
%%%%%%%%%%%%%%%%%%%%%%%%%
[M,I] = max(data,[],1);
[mxv,idx] = max(abs(data(:)));
[r,c,p] = ind2sub(size(data),idx);

r=300;c=300;p=100;
data(r,c,p)
data1_max_loc=[r,c,p,mxv];
% ix=dims(3)/2;
ix=256;
val=1.5e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1e-1;
% a=10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r=68;c=132;p=100;
figure
subplot(1,3,1)
imagesc( squeeze(data(r+1,:,:)).',[-a,a]);
title(strcat('SB,x=',num2str(r)));
xlabel('Y');ylabel('Z');
colorbar

subplot(1,3,2)
imagesc( squeeze(data(:,c+1,:)).',[-a,a]  );    %p
title(strcat('SB,y=',num2str(c)));
xlabel('X');ylabel('Z');
colorbar

subplot(1,3,3)
imagesc( squeeze(data(:,:,p+1)).',[-a,a]);    %p
title(strcat('SB,z=',num2str(p)));
xlabel('X');ylabel('Y');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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