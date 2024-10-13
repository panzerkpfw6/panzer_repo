clear all
close all
addpath('./segymat')
[Data,SegyTraceHeaders,SegyHeader]=ReadSegy('/Volumes/ssd1/Salt_Model_3D/3-D_Salt_Model/ClassicA/salt-a1.segy');
min(Data,[],'all')
max(Data,[],'all')

V=1100;

nx=766;
ny=766;
nz=300;
A=V*ones(nx,ny,nz);
fid=fopen('marmousi_custom.raw','w');
fwrite(fid,A,'float32');
fclose(fid);
clear A

ns=get_ns([nx,ny,nz]);
ss=1
return

function output_size=get_ns(size)
    nx=size(1);
    ny=size(2);
    nz=size(3);
    nx=nx+8;
    ny=ny+8;
    nz=nz+8;
    output_size=nx*ny*nz;
end

