clear all
close all
% fname=['../rtm_munich/data/sismos_16447.raw'];
% fname=['../rtm_munich/data/sismos_91.raw'];
fname=['../stencil-rtm/data/sismos_91.raw'];
fname=['../../data/sismos_20464.raw'];
fname=['../../data/sismos_20475.raw'];

fdir='../../data';
fnames=dir(fullfile(fdir,'*sismos*'));
% fname=fullfile(fdir,fnames(1).name);

nt=1200;nx=128;ny=256;
nt=2200;nx=676;ny=676;
% nt=2200;nx=128;ny=256;
% nt=2200;nx=256;ny=128;

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
imagesc( squeeze(data(310,:,:)).' );
% imagesc( squeeze(data(:,11,:)).' );
clim([-val,val])
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
            data(i,:,k)=A(tmp:(tmp+ny-1));    %x
            tmp=tmp+ny;
        end
    end
end


% function data=read_sismos_orig(fname,dims)
%     ccnt=dims(1)*dims(2)*dims(3);
%     fileID = fopen(fname,'r');
%     A=fread(fileID,ccnt,'float32');
%     fclose(fileID);
%     s = dir(fname);filesize = s.bytes;
% 
%     nx=dims(1);
%     data=nan(dims);
%     tmp=1;
%     for k=1:dims(3)    %z
% %         disp(k)
%         for j=1:dims(2)    %y
%             data(:,j,k)=A(tmp:(tmp+nx-1));    %x
%             tmp=tmp+nx;
%         end
%     end
% end


