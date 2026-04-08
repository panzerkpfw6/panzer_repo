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
%%%%% comparison 1
% root='../../stencil-main';fname= ['SB_1st_abc.raw'];
% root2='../rtm_munich'; fname2=['snapshot_504'];
%%%%% comparison2. SB 1st,SB 2nd.
title1='SB1st';title2='SB2nd';
description_str='diff SB 1st,SB 2nd';
root='../../';  fname= ['snapshot_SB1st_504'];
root2='../../'; fname2=['snapshot_SB2nd_504'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simwave 2nd
description_str='diff SB 2nd,TB 2nd. simwave';title1='SB2nd';title2='TB2nd';
root='../../';  fname= ['snapshot_SB2nd_529'];
root='../../';  fname= ['snapshot_SB2nd_530'];
root2='../../'; fname2=['snapshot_TB2nd_530'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compare "simwave" and "stencil",2nd TB
% description_str='diff TB 2nd. simwave,stencil';title1='simwave';title2='stencil';
% root='../../';  fname= ['snapshot_TB2nd_530'];
% root2='../../../../stencil-main'; fname2=['TB_2nd_abc.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simwave 1st
% description_str='diff SB 1st,TB 1st. simwave';title1='SB1st';title2='TB1st';
% % root='../../';  fname= ['snapshot_SB1st_545'];
% % root2='../../'; fname2=['snapshot_TB1st_546'];
% root='../../';  fname= ['snapshot_SB1st_537'];
% root2='../../'; fname2=['snapshot_TB1st_537'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stencil 2nd
% description_str='diff SB 2nd,TB 2nd.stencil';title1='SB2nd stencil';title2='TB2nd stencil';
% root='../../../../stencil-main';  fname= ['SB_2nd_abc.raw'];
% root2='../../../../stencil-main'; fname2=['TB_2nd_abc.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stencil 1st
% description_str='diff SB 1st,TB 1st.stencil';title1='SB1st stencil';title2='TB1st stencil';
% root='../../../../stencil-main';  fname= ['SB_1st_abc.raw'];
% root2='../../../../stencil-main'; fname2=['TB_1st_abc.raw'];
%%%%%%%%%%%%%%%%%%  compare simwave and stencil. SB 2nd order.!!!!!!!!!!!!!!!!
% description_str='diff SB 2nd,.simwave VS stencil';title1='simwave SB2nd';title2='stencil SB2nd';
% root='../../';  fname= ['snapshot_SB2nd_530'];
% root2='../../../../stencil-main'; fname2=['SB_2nd_abc.raw'];
%%%%%%%%%%%%%  compare simwave and stencil. SB 1st order.!!!!!!!!!!!!!!!!!!!!
% description_str='diff SB 1st,.simwave VS stencil';title1='simwave SB1st';title2='stencil SB1st';
% root='../../';  fname= ['snapshot_SB1st_536'];
% root2='../../../../stencil-main'; fname2=['SB_1st_abc.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% analytical
% description_str='diff SB 2nd,TB 2nd. stencil';title1='SB2nd';title2='Analytical';
% root='../../../../stencil-main'; fname=['SB_2nd_abc.raw'];
% root2='../../../../stencil-main'; fname2=['analytical_sol_529.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stencil-rtm 1st
description_str='diff SB 1st,TB 1st. stencil-rtm';title1='SB1st';title2='TB1st';
root='../../';  fname= ['snapshot_SB1st_505'];
root2='../../'; fname2=['snapshot_TB1st_505'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stencil-rtm 2nd
% description_str='diff SB 2nd,TB 2nd. stencil-rtm';title1='SB2nd';title2='TB2nd';
% root='../../';  fname= ['snapshot_SB2nd_514'];
% root2='../../'; fname2=['snapshot_TB2nd_514'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stencil-rtm SB 1st,2nd
% description_str='diff SB 1st,SB 2nd. stencil-rtm';title1='SB1st';title2='SB2nd';
% root='../../';  fname= ['snapshot_SB1st_519'];
% root2='../../'; fname2=['snapshot_SB2nd_520'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname=fullfile(root,fname);
fname2=fullfile(root2,fname2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add 8 points to each dimension
dims1=[128+8,256+8,512+8];
dims1=[256+8,256+8,256+8];
dims1=[200+8,256+8,160+8];
% dims1=[676+8,676+8,201+8];
ccnt=dims1(1)*dims1(2)*dims1(3);
%%%%%%%%%%%%%%%%%%%%%%%%% check sizes
pwd
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
%%%%%%%%%%%%%%%%%%%%%%%%% data
% data=read_snap(fname,'simwave',dims1);
data=read_snap(fname,'stencil',dims1);
%%%%%%%%%%%%%%%%%%%%%%%%% data2
% data2=read_snap(fname2,'simwave',dims1);
data2=read_snap(fname2,'stencil',dims1);
%%%%%%%%%%%%%%%%%%%%%%%%%
data_diff=data-data2;
1+1e-11-1
C =1+ 0.1.^(+12)-1

min(data_diff,[],'all')
max(data_diff,[],'all')
abs_diff=max(abs(data_diff),[],'all')
RMS_val =rms(data_diff,"all")

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

[M,I] = max(data_diff,[],1);
[mxv,idx] = max(abs(data_diff(:)));
[r2,c2,p2] = ind2sub(size(data_diff),idx);
data_diff_max_loc=[r2,c2,p2,mxv];

% ix=dims(3)/2;
ix=256;
val=1.5e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=5*1e-2;
% a=7;
a=max(abs(data2(:)));
r2=68;c2=132;p2=260;
r2=132;c2=132;p2=260;
r2=132;c2=132;p2=132;
r2=160;c2=140;p2=60;
max(abs(data2(:)))
figure
subplot(1,3,1)
% imagesc( squeeze(data2(r2,:,:)).' );    %p2
imagesc( squeeze(data2(r2,:,:)).' ,[-a,a] );    %p2
title(strcat(title2,',x=',num2str(r2)) );
xlabel('Y');ylabel('Z');
colorbar

subplot(1,3,2)
% imagesc( squeeze(data2(:,c2,:)).' );    %p2
imagesc( squeeze(data2(:,c2,:)).' ,[-a,a] );    %p2
title(strcat(title2,',y=',num2str(c2)));
xlabel('X');ylabel('Z');
colorbar

subplot(1,3,3)
% imagesc( squeeze(data2(:,:,p2)).');    %p2
imagesc( squeeze(data2(:,:,p2)).' ,[-a,a] );    %p2
title(strcat(title2,',z=',num2str(p2)));
xlabel('X');ylabel('Y');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a=max(abs(data(:)));
% r=68;c=132;p=260;
r=r2;c=c2;p=p2;
figure
subplot(1,3,1)
% imagesc( squeeze(data(r+1,:,:)).');
imagesc( squeeze(data(r+1,:,:)).',[-a,a]);
title(strcat(title1,',x=',num2str(r)));
xlabel('Y');ylabel('Z');
colorbar

subplot(1,3,2)
% imagesc( squeeze(data(:,c+1,:)).'  );    %p
imagesc( squeeze(data(:,c+1,:)).',[-a,a]);    %p
title(strcat(title1,',y=',num2str(c)));
xlabel('X');ylabel('Z');
colorbar

subplot(1,3,3)
% imagesc( squeeze(data(:,:,p+1)).'  );    %p
imagesc( squeeze(data(:,:,p+1)).',[-a,a]);    %p
title(strcat(title1,',z=',num2str(p)));
xlabel('X');ylabel('Y');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
title('Difference,\n second order solution')
subplot(1,3,1)
imagesc( squeeze(data_diff(r2,:,:)).' );    %p2
title(strcat(description_str,',x=',num2str(r2)) );
xlabel('Y');ylabel('Z');
colorbar

subplot(1,3,2)
imagesc( squeeze(data_diff(:,c2,:)).' );    %p2
title(strcat(description_str,',y=',num2str(c2)));
xlabel('X');ylabel('Z');
colorbar

subplot(1,3,3)
imagesc(squeeze(data_diff(:,:,p2)) );    %p2
title(strcat(description_str,',z=',num2str(p2)));
xlabel('X');ylabel('Y');
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_diff=max(abs(data_diff),[],'all')
RMS_val =rms(data_diff,"all")
ss=1
%%%%%%%%%%%%%%%%%%   plot trace %%%%%%%%%%%%%%%%%%%
figure;
% plot(data(r2,:,p2));hold on
% plot(data2(r2,:,p2));
% plot(data_diff(r2,:,p2));
% plot(data(:,c2,p2));hold on
% plot(data2(:,c2,p2));
plot(data_diff(:,c2,p2));
legend(title1,title2,'diff');
ss=1
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