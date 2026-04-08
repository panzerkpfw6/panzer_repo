clear all
close all

root='/Volumes/ssd1/snapshots';
% root='~/Downloads';
fname= ['SB_2nd_2530.raw'];

fname2=['SB_2nd_abc_BIG_TEST4_2530.raw'];

fname=fullfile(root,fname);
fname2=fullfile(root,fname2);

% add 8 points to each dimension
dims=[520,520,520];
ccnt=dims(1)*dims(2)*dims(3);

% read file according to dimensions
fileID = fopen(fname,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data = reshape(A,dims);

fileID = fopen(fname2,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data2 = reshape(A,dims);

data=permute(data, [3 2 1]);
data2=permute(data2, [3 2 1]);
close all

ix=260;
figure
imagesc( squeeze(data(:,:,ix)).'  );
val=1.5;
clim([-val,val]);title(fname);
colorbar

figure
imagesc( squeeze(data2(:,:,ix)).' );
clim([-val,val]);
title(fname2);
colorbar

diff_slice2=data(:,:,ix)-data2(:,:,ix);

figure
imagesc(diff_slice2);
% clim([-val,val]);
title('diff 2');
colorbar

data_diff=data-data2;
figure
imagesc( squeeze(data_diff(:,:,ix)).' );
clim([-val,val]);
title(fname2);
colorbar

for ix=1:100:500
    figure
    imagesc( squeeze(data_diff(:,:,ix)).'  );
    val=1.5;
    clim([-val,val]);title(fname);
    colorbar
end

min(diff_slice2,[],'all')
max(diff_slice2,[],'all')
min(data2,[],'all')
max(data2,[],'all')

figure
imagesc( squeeze(data_diff(:,:,260)).'  );
val=1.5;
clim([-val,val]);title(fname);
colorbar

figure
imagesc( squeeze(data_diff(:,260,:)).'  );
val=1.5;
clim([-val,val]);title(fname);
colorbar

figure
imagesc( squeeze(data_diff(260,:,:)).'  );
val=1.5;
clim([-val,val]);title(fname);
colorbar

matrix=squeeze(data_diff(:,:,260)).';
profile=matrix(:,113);

figure;
plot(profile);

figure;
plot(matrix(:,250));

ss=1