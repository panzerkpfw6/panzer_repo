clear all
close all

root='/Volumes/ssd1/snapshots/check_abc';
fname= ['SB_1st_abc.raw'];
fname2=['SB_2nd_abc.raw'];
fname3=['TB_1st_abc.raw'];
fname4=['TB_2nd_abc.raw'];

fname=fullfile(root,fname);
fname2=fullfile(root,fname2);
fname3=fullfile(root,fname3);
fname4=fullfile(root,fname4);

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

fileID = fopen(fname3,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data3 = reshape(A,dims);

fileID = fopen(fname4,'r');
A=fread(fileID,ccnt,'float32');
fclose(fileID);
data4 = reshape(A,dims);

data=permute(data, [3 2 1]);
data2=permute(data2, [3 2 1]);
data3=permute(data3, [3 2 1]);
data4=permute(data4, [3 2 1]);


close all

ix=260;
% plotting
figure
imagesc( squeeze(data(:,:,ix)).'  );
val=1.5;
val=2;
clim([-val,val]);title(fname);
colorbar

figure
imagesc( squeeze(data2(:,:,ix)).' );
clim([-val,val]);
title(fname2);
colorbar

figure
imagesc( squeeze(data3(:,:,ix)).' );
clim([-val,val]);title(fname3);
colorbar

figure
imagesc( squeeze(data4(:,:,ix)).' );
clim([-val,val]);
title(fname4);
colorbar

ss=1

diff_slice2=data(:,:,ix)-data3(:,:,ix);
figure
imagesc(diff_slice2);
clim([-val,val]);
title('difference SB 1st, TB 1st ');
colorbar

diff_slice3=data(:,:,ix)-data2(:,:,ix);
figure
imagesc(diff_slice3);
clim([-val,val]);
title('difference SB 1st, SB 2nd ');
colorbar

diff_slice4=data3(:,:,ix)-data4(:,:,ix);
figure
imagesc(diff_slice4);
clim([-val,val]);
title('difference TB 1st, TB 2nd ');
colorbar

diff_slice5=data2(:,:,ix)-data4(:,:,ix);
figure
imagesc(diff_slice5);
clim([-val,val]);
title('difference SB 2nd, TB 2nd ');
colorbar
ss=1

% 
% data_diff=data-data2;
% figure
% imagesc( squeeze(data_diff(:,:,ix)).' );
% clim([-val,val]);
% title(fname2);
% colorbar
% 
% for ix=1:100:500
%     figure
%     imagesc( squeeze(data_diff(:,:,ix)).'  );
%     val=1.5;
%     clim([-val,val]);title(fname);
%     colorbar
% end
% 
% min(diff_slice2,[],'all')
% max(diff_slice2,[],'all')
% min(data2,[],'all')
% max(data2,[],'all')
% 
% figure
% imagesc( squeeze(data_diff(:,:,260)).'  );
% val=1.5;
% clim([-val,val]);title(fname);
% colorbar
% 
% figure
% imagesc( squeeze(data_diff(:,260,:)).'  );
% val=1.5;
% clim([-val,val]);title(fname);
% colorbar
% 
% figure
% imagesc( squeeze(data_diff(260,:,:)).'  );
% val=1.5;
% clim([-val,val]);title(fname);
% colorbar

ss=1