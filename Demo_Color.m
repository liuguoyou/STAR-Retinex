% Retinex theory is a color perception model of human vision and is
% used to remove illumination effects in images. The primary goal
% of Retinex is to decompose the observed images into illumination
% and reflectance.

clc;clear;
%%% methods
addpath(genpath('methods'));
%%% test images
ext         =  {'*.tif','*.jpg','*.jpeg','*.JPG','*.png','*.bmp'};
im_dir   =  [];
for i = 1 : length(ext)
    im_dir = cat(1,im_dir, dir(fullfile('Color', ext{i})));
end
im_num = length(im_dir);

i=2;
name = regexp(im_dir(i).name, '\.', 'split');
Im=im2double( imread(fullfile('Color', im_dir(i).name)) );

% JIEP ICCV2017
method = 'JieP';
ccIm = zeros(size(Im));
for ch = 1:size(Im,3)
    S = Im(:,:,ch);
    [I, R] = jiep(S);
    ccIm(:,:,ch) = R;
end
imwrite(ccIm, [name{1} '_' method '.png']);

% STAR
method = 'STAR';
alpha = 0.001;
beta = 0.0001;
ccIm = zeros(size(Im));
for pI = [1.5]
    for pR = [0.5]
        for ch = 1:size(Im,3)
            S = Im(:,:,ch);
            [I, R] = STAR(Im, alpha, beta, pI, pR);
            ccIm(:,:,ch) = R;
        end
        imwrite(ccIm(hsv), [name{1} '_' method '.png']);
    end
end


% WVM CVPR2016
method = 'WVM';
Im = 255*Im;
ccIm = zeros(size(Im));
for ch = 1:size(Im,3)
    S = Im(:,:,ch);
    c_1 = 0.01; c_2 = 0.1; lambda = 1;     % set parameters
    epsilon_stop = 1e-3;  % stopping criteria
    [ R, I, epsilon_R, epsilon_L ] = WVM_CVPR2016( S, c_1, c_2, lambda, epsilon_stop );
    ccIm(:,:,ch) = R;
end
imwrite(ccIm, [name{1} '_' method '.png']);