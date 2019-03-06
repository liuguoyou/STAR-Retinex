im = '1';
h  = 160;
w  = 80;
s  = 35;
f  = 3;
lr = 1;
scale = 1/3;
ii = 1;

%% input
image = imread(['/home/csjunxu/Paper/Enhancement/Dataset/Images_LowLight/' im '.bmp']);
[hh, ww, cc] = size(image);
ll = min(hh,ww);
image = imresize(image(1:ll,ii:ll+ii-1,:), scale);
[ outputimage ] = boxandresize( image, h,w,s, f,lr);
imname = sprintf('%s_%s.png','rs',im);
imwrite(outputimage,imname,'png');

%% NPE TIP2013
image = imread(['/home/csjunxu/Paper/Enhancement/Results_LowLight/NPE_TIP2013/NPE_TIP2013_' im '.jpg']);
image = imresize(image(1:ll,ii:ll+ii-1,:), scale);
[ outputimage ] = boxandresize( image, h,w,s, f,lr);
imname = sprintf('%s_%s_NPE.png','rs',im);
imwrite(outputimage,imname,'png');

%% WVM CVPR2016
% ours
image = imread(['/home/csjunxu/Paper/Enhancement/Results_LowLight/WVM_CVPR2016/WVM_CVPR2016_' im '.jpg']);
image = imresize(image(1:ll,ii:ll+ii-1,:), scale);
[ outputimage ] = boxandresize( image, h,w,s, f,lr);
imname = sprintf('%s_%s_WVM.png','rs',im);
imwrite(outputimage,imname,'png');

%% JieP ICCV 2017
image = imread(['/home/csjunxu/Paper/Enhancement/Results_LowLight/JieP_ICCV2017/JieP_ICCV2017_' im '.jpg']);
image = imresize(image(1:ll,ii:ll+ii-1,:), scale);
[ outputimage ] = boxandresize( image, h,w,s, f,lr);
imname = sprintf('%s_%s_JieP.png','rs',im);
imwrite(outputimage,imname,'png');

%% LIME TIP2017
% ours
image = imread(['/home/csjunxu/Paper/Enhancement/Results_LowLight/LIME_TIP2017/LIME_TIP2017_' im '.jpg']);
image = imresize(image(1:ll,ii:ll+ii-1,:), scale);
[ outputimage ] = boxandresize( image, h,w,s, f,lr);
imname = sprintf('%s_%s_LIME.png','rs',im);
imwrite(outputimage,imname,'png');

%% STAR
image = imread(['/home/csjunxu/Paper/Enhancement/Results_LowLight/STAR/STAR_10_aIpI=2_RpR=1_alpha=0.001_beta=0.0001_' im '.png']);
image = imresize(image(1:ll,ii:ll+ii-1,:), scale);
[ outputimage ] = boxandresize( image, h,w,s, f,lr);
imname = sprintf('%s_%s_STAR.png','rs',im);
imwrite(outputimage,imname,'png');


