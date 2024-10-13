clear all
close all

root='/Volumes/ssd1/snapshots';
fname=['SB_1st_2301_512_751_1500.raw'];
fname=fullfile(root,fname);
dims=[2301,512,751];

% add 8 points to each dimension
dims=dims+8;
ccnt=dims(1)*dims(2)*dims(3);

% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data= reshape(A,[dims(3),dims(2),dims(1)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; imagesc( squeeze(data(354,:,:)) );colorbar
figure; imagesc( squeeze(data(:,304,:)) );colorbar
figure; imagesc( squeeze(data(:,:,254)) );colorbar

[M,I]=max(data,[],'all');
[i1,i2,i3] = ind2sub(size(data),I);

figure;imagesc( squeeze(data(:,i2,:)).' );colorbar

figure;imagesc( squeeze(data(i1,:,:)).' );colorbar

figure;imagesc( squeeze(data(:,:,i3)).' );colorbar
