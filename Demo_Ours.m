clc;
clear;
Original_image_dir  =    '/home/csjunxu/Paper/Enhancement/Dataset/LowLightImages/'; % put clean images in this folder
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.bmp');
im_dir  = dir(fpath);
im_num = length(im_dir);

% metrics
addpath('metrics');

% write_dir  = ['/home/csjunxu/Github/data/Ultrasound/'];
% write_dir  = ['/home/csjunxu/Github/Segmentation-master/'];
write_dir  = ['/home/csjunxu/Paper/Enhancement/Results_LowLight/'];
if ~isdir(write_dir)
    mkdir(write_dir);
end

gamma=2.2;
alpha = 0.001;
beta = 0.0001;
for pI = [.1:.1:2]
    for pR = [.3:.1:2]
        NIQEs = zeros(im_num,1);
        LOEs = zeros(im_num,1);
        for i = 1:im_num
            Im=im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            [I, R] = enhancer(Im, alpha, beta, pI, pR);
            hsv = rgb2hsv(Im);
            I_gamma = I.^(1/gamma);
            S_gamma = R .* I_gamma;
            hsv(:,:,3) = S_gamma;
            eIm = hsv2rgb(hsv);
            NIQEs(i) = niqe(eIm);
            LOEs(i) = LOE(eIm, Im);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
            % imwrite(enhance, [write_dir im_dir(i).name '_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
            %    num2str(alpha) '_beta=' num2str(beta) '.jpg'])
        end
        matname = [write_dir '/Our_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
            num2str(alpha) '_beta=' num2str(beta) '.mat'];
        mNIQEs = mean(NIQEs);
        mLOEs = mean(LOEs);
        fprintf('Mean NIQE = %2.4f, Mean LOE = %2.4f\n', mNIQEs, mLOEs);
        save(matname, 'NIQEs', 'mNIQEs', 'LOEs', 'mLOEs');
    end
end
% subplot(1,2,1); imshow(im); title('Input');
% subplot(1,2,2); imshow(enhance); title('Illumination Adjustment');