clc;clear;
% test images
Original_image_dir  =    '/home/csjunxu/Paper/Enhancement/Dataset/LowLightImages/'; % put clean images in this folder
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.bmp');
im_dir  = dir(fpath);
im_num = length(im_dir);

% metrics
addpath('metrics');
metrics = {'LOE', 'NIQE'};
% methods
methods = {'JieP_ICCV2017', 'HE'};

for m = 1%:length(methods)
    method = methods{m};
    % write_mat_dir  = ['/home/csjunxu/Github/data/Ultrasound/'];
    % write_mat_dir  = ['/home/csjunxu/Github/Segmentation-master/'];
    write_mat_dir = '/home/csjunxu/Paper/Enhancement/Results_LowLight/';
    write_img_dir = [write_mat_dir method '/'];
    if ~isdir(write_img_dir)
        mkdir(write_img_dir);
    end
    NIQEs = zeros(im_num,1);
    LOEs = zeros(im_num,1);
    for i = 1:im_num
        S = regexp(im_dir(i).name, '\.', 'split');
        Im=im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
        if strcmp(method, 'JieP_ICCV2017') == 1
            gamma=2.2;
            [I, R] = jiep(Im);
            hsv = rgb2hsv(Im);
            I_gamma = I.^(1/gamma);
            S_gamma = R .* I_gamma;
            hsv(:,:,3) = S_gamma;
            eIm = hsv2rgb(hsv);
        elseif strcmp(method, 'HE') == 1
            [eIm, ~] = histeq(Im);
        end
        NIQEs(i) = niqe(eIm);
        LOEs(i) = LOE(eIm, Im);
        fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        imwrite(eIm, [write_img_dir method '_' S{1} '.jpg']);
    end
    matname = [write_mat_dir method '.mat'];
    mNIQEs = mean(NIQEs);
    mLOEs = mean(LOEs);
    fprintf('Mean NIQE = %2.4f, Mean LOE = %2.4f\n', mNIQEs, mLOEs);
    save(matname, 'NIQEs', 'mNIQEs', 'LOEs', 'mLOEs');
end
% subplot(1,2,1); imshow(im); title('Input');
% subplot(1,2,2); imshow(enhance); title('Illumination Adjustment');