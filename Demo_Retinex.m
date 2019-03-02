% Retinex theory is a color perception model of human vision and is
% used to remove illumination effects in images. The primary goal
% of Retinex is to decompose the observed images into illumination
% and reflectance.

clc;clear;
%%% choose test dataset
datasets = {'LowLight', 'NPE', 'VV', 'NASA', 'LDR'};
Testset = datasets{1}; % select test dataset
Test_dir  = fullfile('/home/csjunxu/Paper/Enhancement/Dataset', ['Images_' Testset]);
%%% read images
ext         =  {'*.jpg','*.jpeg','*.JPG','*.png','*.bmp'};
im_dir   =  [];
for i = 1 : length(ext)
    im_dir = cat(1,im_dir, dir(fullfile(Test_dir,ext{i})));
end
im_num = max(8,length(im_dir));

name = regexp(im_dir(im_num).name, '\.', 'split');
Im=im2double( imread(fullfile(Test_dir, im_dir(im_num).name)) );
write_dir = '/home/csjunxu/Paper/Enhancement/Results_Retinex/';
if ~isdir(write_dir)
    mkdir(write_dir);
end
imwrite(Im, [write_dir name{1} '.png'])
% STAR
method = 'STAR';
alpha = 0.001;
beta = 0.0001;
for pI = [1.5:.1:2]
    for pR = [0.5:-.1:.1]
        [I, R] = STAR(Im, alpha, beta, pI, pR);
        hsv = rgb2hsv(Im);
        subplot(2,2,1); imshow(I);  title('Illumination (Gray)');
        imwrite(I, [write_dir name{1} '_I_Gray_' method '_pI' num2str(pI) '_pR' num2str(pR) '.png'])
        hsv(:,:,3) = I;
        subplot(2,2,2); imshow(hsv2rgb(hsv));  title('Illumination (RGB)');
        imwrite(hsv2rgb(hsv), [write_dir name{1} '_I_RGB_' method '_pI' num2str(pI) '_pR' num2str(pR) '.png'])
        subplot(2,2,3); imshow(I);  title('Reflectance (Gray)');
        imwrite(R, [write_dir name{1} '_R_Gray_' method '_pI' num2str(pI) '_pR' num2str(pR) '.png'])
        hsv(:,:,3) = R;
        subplot(2,2,4); imshow(hsv2rgb(hsv));  title('Reflectance (RGB)');
        imwrite(hsv2rgb(hsv), [write_dir name{1} '_R_RGB_' method '_pI' num2str(pI) '_pR' num2str(pR) '.png'])
    end
end

% JIEP ICCV2017
method = 'JieP';
[I, R] = jiep(Im);
hsv = rgb2hsv(Im);
subplot(2,2,1); imshow(I);  title('Illumination (Gray)');
imwrite(I, [write_dir name{1} '_I_Gray_' method '.png'])
hsv(:,:,3) = I;
subplot(2,2,2); imshow(hsv2rgb(hsv));  title('Illumination (RGB)');
imwrite(hsv2rgb(hsv), [write_dir name{1} '_I_RGB_' method '.png'])
subplot(2,2,3); imshow(I);  title('Reflectance (Gray)');
imwrite(R, [write_dir name{1} '_R_Gray_' method '.png'])
hsv(:,:,3) = R;
subplot(2,2,4); imshow(hsv2rgb(hsv));  title('Reflectance (RGB)');
imwrite(hsv2rgb(hsv), [write_dir name{1} '_R_RGB_' method '.png'])
%subplot(2,2,1); imshow(Im);  title('Input');