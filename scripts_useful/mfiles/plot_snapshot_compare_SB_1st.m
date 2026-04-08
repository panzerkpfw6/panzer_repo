clear all
close all

root='/Volumes/ssd1/snapshots';
root='../../stencil-sismos/SB_abc';

root='../../stencil-ABCs-left_side_fix/raws';
root='../../tb_abc_fix';
% root='~/Downloads';
fname= ['SB_1st.raw'];

% fname2=['SB_1st_2.raw'];
% fname2=['SB_2nd_abc_BIG_TEST4_2530.raw'];
fname2=['SB_1st_abc.raw'];

fname= ['TB_1st_2530.raw'];
fname2=['TB_1st_abc.raw'];

fname= ['SB_1st_abc.raw'];
fname2=['SB_2nd_abc.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% root='../../stencil-tb_abc_fix/TB_abc';
% % root='../../stencil-tb_abc_fix/raws';
% root='../../stencil-tb_abc_fix/TB_abc';
% fname= ['TB_1st_abc.raw'];
% fname2=['TB_1st_abc.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='../../main_test_branch/raws';root2=root;
fname= ['SB_1st.raw'];
fname2=['TB_1st.raw'];
% fname2=['TB_1st_abc.raw'];
% fname2=['TB_2nd.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% root='../../main_24_09_23/raws';root2=root;
% fname= ['SB_1st.raw'];
% fname2=['TB_1st.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% root='../../main_test_branch/raws';fname= ['SB_1st.raw'];
% root2='../../main_24_09_23/raws';fname2=['SB_1st.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='../../VALIDITY_TEST';root2=root;
fname= ['SB_1st.raw'];
fname2=['TB_1st.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% root='../../MAIN2/raws';root2=root;
% fname= ['SB_1st.raw'];
% fname2=['TB_1st.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='../../tb-sismos/raws';root2=root;
fname= ['SB_1st_abc.raw'];
fname2=['TB_1st_abc.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='../../stencil-tb_abc_fix_dead_branch';root2=root;
fname= ['SB_1st_abc.raw'];
fname2=['TB_1st_abc.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='../../tb_abc_fix';root2=root;
fname= ['SB_1st_abc.raw'];
fname2=['TB_1st_abc.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root='../../tb_abc_fix';root2=root;
fname= ['SB_2nd.raw'];
fname2=['TB_2nd.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% root ='../../tb-sismos/SB_abc';
% fname= ['SB_1st_abc.raw'];
% root2='../../tb-sismos/TB_abc';
% fname2=['TB_1st_abc.raw'];
% % root2='../../tb-sismos/TB';
% % fname2=['TB_1st.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root ='../../tb-sismos';
fname= ['SB_1st.raw'];
root2='../../tb-sismos';
fname2=['TB_1st.raw'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname=fullfile(root,fname);
fname2=fullfile(root2,fname2);

% add 8 points to each dimension
dims=[520,520,520];
% dims=[264,264,264];
% dims=[264,264,264];
% dims=[308,308,308];
% dims=[2056,2056,520];
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

data=permute(data,[3 2 1]);
data2=permute(data2,[3 2 1]);
data3=abs(data2);

data_diff=data-data2;
min(data_diff,[],'all')
max(data_diff,[],'all')

min(data,[],'all')
max(data,[],'all')

min(data2,[],'all')
max(data2,[],'all')

close all

min(data,[],'all')
max(data,[],'all')

max(data);
[M,I] = max(data,[],1);
[mxv,idx] = max(data(:));
[r,c,p] = ind2sub(size(data),idx);
data(r,c,p)

[M,I] = max(data2,[],1);
[mxv,idx] = max(data2(:));
[r2,c2,p2] = ind2sub(size(data2),idx);


max(abs(data_diff(:)))
[M,I] = max(data_diff,[],1);
[mxv,idx] = max(data_diff(:));
[r,c,p] = ind2sub(size(data_diff),idx);
data_diff(r,c,p)

figure
imagesc( squeeze(data_diff(:,:,ix)).' );
title('diff 2nd SB TB');
colorbar

ix=dims(3)/2;
% ix=133;
% ix=10;
val=1.5e-4;
% val=1.0;

% ix=r;
figure
imagesc( squeeze(data2(:,:,ix)).'  );
% caxis([-val,val]);
title(fname2);
colorbar

figure
imagesc( squeeze(data(:,:,ix)).'  );
% caxis([-val,val]);
title(fname);
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
imagesc( squeeze(data2(:,:,10)).'  );
% caxis([-val,val]);
title(fname2);
colorbar

figure
imagesc( squeeze(data(:,:,101)).'  );
% caxis([-val,val]);
title(fname);
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
% imagesc( squeeze(data(:,:,ix)).'  );
imagesc( squeeze(data(:,ix,:)).'  );
caxis([-val,val]);
title(fname);
colorbar

figure
imagesc( squeeze(data(ix,:,:)).'  );
caxis([-val,val]);
title(fname);
colorbar

figure
imagesc( squeeze(data(:,:,ix)).' );
caxis([-val,val]);
title(fname);
colorbar

figure
imagesc( squeeze(data2(:,:,ix)).' );
caxis([-val,val]);
title(fname2);
colorbar

figure
imagesc( squeeze(data_diff(:,:,ix)).' );
% imagesc( squeeze(data_diff(:,ix,:)).' );
% imagesc( squeeze(data_diff(ix,:,:)).' );
caxis([-val,val]);
title('diff sb 1st and 2nd');
colorbar

figure
imagesc( squeeze(data_diff(:,:,ix)).' );
title('diff 2nd SB TB');
colorbar

figure
plot( squeeze(data(:,ix,ix)).' );hold on
plot( squeeze(data2(:,ix,ix)).' );
% plot( squeeze(data_diff(:,ix,ix)).' );
title('sb tb 1st');
legend('SB','MWD-TB')
colorbar

figure
% plot( squeeze(data(:,ix,ix)).' );hold on
% plot( squeeze(data2(:,ix,ix)).' );
plot( squeeze(data_diff(:,ix,ix)).' );
title('diff sb tb 1st');
legend('diff')
colorbar

figure
% plot( squeeze(data(:,ix,ix)).' );hold on
% plot( squeeze(data2(:,ix,ix)).' );
plot( squeeze(data_diff(:,ix,ix)).' );
title('diff sb 1st and tb 1st');
legend
colorbar

min(data,[],'all')
max(data,[],'all')

for ix=1:100:500
    figure
    imagesc( squeeze(data_diff(:,:,ix)).'  );
    val=1.5;
    caxis([-val,val]);title(fname);
    colorbar
end

min(diff_slice2,[],'all')
max(diff_slice2,[],'all')
min(data2,[],'all')
max(data2,[],'all')

figure
imagesc( squeeze(data_diff(:,:,dims(3)/2)).'  );
val=1.5;
caxis([-val,val]);title(fname);
colorbar

figure
imagesc( squeeze(data_diff(:,dims(2)/2,:)).'  );
val=1.5;
caxis([-val,val]);title(fname);
colorbar

figure
imagesc( squeeze(data_diff(dims(1)/2,:,:)).'  );
val=1.5;
caxis([-val,val]);title(fname);
colorbar

ss=1