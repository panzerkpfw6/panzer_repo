clear all
close all
% code to compare source implementation in different codes
times=[];
%%  stencil-rtm sources
fname=['../../data/firstsrc.raw'];
fname2=['../../data/secondsrc.raw'];

% some file size calculations
s = dir(fname);
s2 = dir(fname2);
nt1= s.bytes/4;
nt2= s2.bytes/4;
times=[times,nt1,nt2];

% read file according to dimensions
fileID = fopen(fname,'r');
simwave_1st=fread(fileID,nt1,'float');
fclose(fileID);

fileID = fopen(fname2,'r');
simwave_2nd=fread(fileID,nt2,'float');
fclose(fileID);
%%  stencil sources
fname= ['../../../../stencil-main/data/stencil_src_1st.raw'];
fname2=['../../../../stencil-main/data/stencil_src_2nd.raw'];

% some file size calculations
s = dir(fname);
s2 = dir(fname2);
nt1= s.bytes/4;
nt2= s2.bytes/4;
times=[times,nt1,nt2];

% read file according to dimensions
fileID = fopen(fname,'r');
stencil_1st=fread(fileID,nt1,'float');
fclose(fileID);

fileID = fopen(fname2,'r');
stencil_2nd=fread(fileID,nt2,'float');
fclose(fileID);
%%
nt=min(times);
t=(1:nt)-1;

figure 
hold on;
plot(t,simwave_1st(1:nt),'Color','r','LineWidth',1);
plot(t,simwave_2nd(1:nt),'Color','r','LineWidth',2);
plot(t,stencil_1st(1:nt),'--','Color','b','LineWidth',1);
plot(t,stencil_2nd(1:nt),'--','Color','b','LineWidth',2);
legend('simwave_1st','simwave_2nd','stencil_1st','stencil_2nd');
%%
figure 
hold on;
plot(t,simwave_1st(1:nt),'Color','r','LineWidth',1);
plot(t,stencil_1st(1:nt),'--','Color','b','LineWidth',1);
legend('simwave_1st','stencil_1st');
%%
figure 
hold on;
plot(t,simwave_2nd(1:nt),'Color','r','LineWidth',2);
plot(t,stencil_2nd(1:nt),'--','Color','b','LineWidth',2);
legend('simwave_2nd','stencil_2nd');

ss=1;