clear all
close all
V=1500;
% record file according to dimensions
ns=get_v(512);

A=1500*ones(512,512,512);
fid=fopen('benchmark_512_3d.raw','w');
fwrite(fid,A,'float32');
fclose(fid);
clear A

% B=1500*ones(1024,1024,512);
% fid=fopen('benchmark_1024_3d.raw','w');
% fwrite(fid,B,'float32');
% fclose(fid);
% clear B
% 
% C=1500*ones(2048,2048,512);
% fid=fopen('benchmark_2048_3d.raw','w');
% fwrite(fid,C,'float32');
% fclose(fid);
% clear C
% ss=1
return

% save_step=10000
% ns2=rem(ns,save_step);
% ns1=ns-ns2
% V=1500*ones(save_step);
% fid=fopen('benchmark_512_512_512.raw','w');
% for i=1:floor(ns/save_step)
% % for i=1:5
%     disp(i)
%     fwrite(fid,V,'float32');
% end
% V=1500*ones(ns2);
% fwrite(fid,V,'float32');
% fclose(fid);
% return
% 
% nx=512
% nx=nx+8;
% ny=nx;
% nz=nx;
% % V=1500*ones(nx*ny*nz)
% ss=1

% fid=fopen('benchmark_1024_1024_1024.raw','w');
% fwrite(fid,V,'float32');
% fclose(fid);

V=1500*ones(1000);
fid=fopen('benchmark_512_512_512.raw','w');
for i=1:100
    fwrite(fid,V,'float32');
end
fclose(fid);

V=get_v(1024);
fid=fopen('benchmark_1024_1024_1024.raw','w');
fwrite(fid,V,'float32');
fclose(fid);

V=get_v(2056);
fid=fopen('benchmark_2056_2056_2056.raw','w');
fwrite(fid,V,'float32');
fclose(fid);

function output_size=get_v(nx)
% function V=get_v(nx)
    nx=nx+8;
    ny=nx;
    nz=512;
%     V=1500*ones(nx*ny*nz);
    output_size=nx*ny*nz;
end

