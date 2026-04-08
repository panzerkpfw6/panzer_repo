clear all
close all
% fname=['../rtm_munich/data/sismos_16447.raw'];
% fname=['../rtm_munich/data/sismos_91.raw'];
fname=['../stencil-rtm/data/sismos_91.raw'];
fname=['../../data/sismos_sb/sismos_1.raw'];
fname=['../../data/sismos_sb/sismos_41100.raw'];
fname=['../../data/sismos_61.raw'];

% fname=['../../data/filtered_real_data/sismos_61.raw'];
nt=1200;nx=128;ny=256;
nt=200;nx=676;ny=676;

nt=2200;nx=128;ny=256;

nt=505;nx=256;ny=256;

nt=4000;nx=676;ny=676;

dims=[nx,ny,nt];
rcv_len=nx*ny;
ccnt=rcv_len*nt;
%%
% some file size calculations
s = dir(fname);         
filesize = s.bytes 
domain_size=ccnt*4;
format longG
disp(strcat('theoretical file size= ',int2str(domain_size)))
disp(strcat('true file size= ',int2str(filesize)))
%%
data=read_sismos(fname,dims);
[M,I] = max(data,[],1);
[mxv,idx] = max(abs(data(:)));
[r,c,p] = ind2sub(size(data),idx);
data(r,c,p)
data1_max_loc=[r,c,p,mxv];

%%  % plotting
% figure
% % imagesc( squeeze(data(:,:,11)).' );
% val=1e-12;
% % clim([-val,val])
% colorbar
val=1e-1;
figure
% imagesc( squeeze(data(:,123,:)).',[-val,val] );
imagesc( squeeze(data(:,136,:)).' );
% clim([-val,val])
colorbar

min(data,[],'all')
max(data,[],'all')

function data=read_sismos(fname,dims)
    ccnt=dims(1)*dims(2)*dims(3);
    fileID = fopen(fname,'r');
    A=fread(fileID,ccnt,'float32');
    fclose(fileID);
    s = dir(fname);filesize = s.bytes;
    
    nx=dims(1);
    ny=dims(2);
    nt=dims(3);
    data=nan(dims);
    tmp=1;
    for k=1:nt    %z
%         disp(k)
        for i=1:nx    %x
            data(i,:,k)=A(tmp:(tmp+ny-1));    %y
            tmp=tmp+ny;
        end
    end
end
