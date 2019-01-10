clc;
clear;
Original_image_dir  =    '/home/csjunxu/Paper/Enhancement/Dataset/LowLightImages/'; % put clean images in this folder
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.bmp');
im_dir  = dir(fpath);
im_num = length(im_dir);

% metrics
addpath('metrics');
metrics = {'LOE', 'NIQE'};

method = 'Our';

% write_mat_dir  = ['/home/csjunxu/Github/data/Ultrasound/'];
% write_mat_dir  = ['/home/csjunxu/Github/Segmentation-master/'];
write_mat_dir = '/home/csjunxu/Paper/Enhancement/Results_LowLight/';
write_img_dir = [write_mat_dir method '/'];
if ~isdir(write_img_dir)
    mkdir(write_img_dir);
end

gamma=2.2;
alpha = 0.001;
beta = 0.0001;
for pI = [.1:.1:2]
    for pR = [.1:.1:2]
        NIQEs = zeros(im_num,1);
        LOEs = zeros(im_num,1);
        for i = 1:im_num
            name = regexp(im_dir(i).name, '\.', 'split');
            Im=im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            [I, R] = enhancer(Im, alpha, beta, pI, pR);
            hsv = rgb2hsv(Im);
            I_gamma = I.^(1/gamma);
            S_gamma = R .* I_gamma;
            hsv(:,:,3) = S_gamma;
            eIm = hsv2rgb(hsv);
            NIQEs(i) = niqe(uint8(eIm*255)); %range: 0~255 -> uint8
            LOEs(i) = LOE(eIm, Im); % range: 0~1
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
            % imwrite(enhance, [write_img_dir method '_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
            %    num2str(alpha) '_beta=' num2str(beta) '_' name{1} '.jpg'])
        end
        matname = [write_mat_dir '/Our_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
            num2str(alpha) '_beta=' num2str(beta) '.mat'];
        mNIQEs = mean(NIQEs);
        mLOEs = mean(LOEs);
        fprintf('Mean NIQE = %2.4f, Mean LOE = %2.4f\n', mNIQEs, mLOEs);
        save(matname, 'NIQEs', 'mNIQEs', 'LOEs', 'mLOEs');
    end
end
% subplot(1,2,1); imshow(im); title('Input');
% subplot(1,2,2); imshow(enhance); title('Illumination Adjustment');