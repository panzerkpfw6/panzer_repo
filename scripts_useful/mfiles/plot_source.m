clear all
close all
nt=504;
fname=['../rtm_munich/data/src.raw'];

% some file size calculations
s = dir(fname);         
filesize = s.bytes ;

% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,nt,'float');
fclose(fileID);

figure 
plot(A)