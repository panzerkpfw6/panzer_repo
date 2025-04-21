clear all
close all
% fname=['../rtm_munich/data/sismos_16447.raw'];
% fname=['../rtm_munich/data/sismos_91.raw'];
fname=['../stencil-rtm/data/sismos_91.raw'];
fname=['../../data/sismos_1.raw'];
fname=['../../data/sismos_20958.raw'];

fdir='../../data';
fdir='../../data/filtered_real_data';
fnames=dir(fullfile(fdir,'*sismos*'));
fname=fullfile(fdir,fnames(1).name);

nt=2200;nx=676;ny=676;

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
val=1e-4;
figure
% imagesc( squeeze(data(:,123,:)).',[-val,val] );
% imagesc( squeeze(data(128,:,:)).' );
imagesc( squeeze(data(:,310,:)).' );
clim([-val,val])
colorbar

min(data,[],'all')
max(data,[],'all')
ss=1

function data=read_sismos(fname,dims)
    ccnt=dims(1)*dims(2)*dims(3);
    fileID = fopen(fname,'r');
    A=fread(fileID,ccnt,'float32');
    fclose(fileID);
    s = dir(fname);filesize = s.bytes;
    
    nx=dims(1);
    data=nan(dims);
    tmp=1;
    for k=1:dims(3)    %z
%         disp(k)
        for j=1:dims(2)    %y
            data(:,j,k)=A(tmp:(tmp+nx-1));    %x
            tmp=tmp+nx;
        end
    end
end


