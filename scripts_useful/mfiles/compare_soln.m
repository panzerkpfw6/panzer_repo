clear all
close all
dims=[128+8,256+8,512+8];
n=dims(1)*dims(2)*dims(3);
fname1=['../rtm_munich/solution_before_2'];
fname2=['../rtm_munich/solution_p_sveep1'];

fname1=['../rtm_munich/solution_before_2'];
fname2=['../rtm_munich/solution_p_sveep1'];

fname1=['../rtm_munich/solution_before_150'];
fname2=['../rtm_munich/solution_p_sveep150'];

% some file size calculations
s = dir(fname1);         
filesize = s.bytes ;

% read file according to dimensions
fileID = fopen(fname1,'r');
A1=fread(fileID,n,'float');
fclose(fileID);

fileID = fopen(fname2,'r');
A2=fread(fileID,n,'float');
fclose(fileID);

min(abs(A1),[],'all')
min(abs(A2),[],'all')

min(A1)
max(A1)
min(A2)
max(A2)

diff=A1-A2;
min(diff)
max(diff)

figure 
plot(A1-A2)

figure 
plot(A1)

