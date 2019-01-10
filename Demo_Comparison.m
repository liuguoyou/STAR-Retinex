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
methods = {'Li_TIP2018', 'JieP_ICCV2017', 'WVM_CVPR2016', 'HE'};
addpath('methods');

for m = 1:length(methods)
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
        name = regexp(im_dir(i).name, '\.', 'split');
        if  strcmp(method, 'Li_TIP2018') == 1
            Im = double(imread(fullfile(Original_image_dir, im_dir(i).name)));
            
            para.epsilon_stop_L = 1e-3;
            para.epsilon_stop_R = 1e-3;
            para.epsilon = 10/255;
            para.u = 1;
            para.ro = 1.5;
            para.lambda = 5;
            para.beta = 0.01;
            para.omega = 0.01;
            para.delta = 10;
            
            gamma = 2.2;
            [R, L, N] = Li_TIP2018(Im, para);
            eIm = R.*L.^(1/gamma);
            imwrite(eIm/255, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(eIm/255);
            LOEs(i) = LOE(eIm/255, Im/255);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        elseif strcmp(method, 'JieP_ICCV2017') == 1
            Im=im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            gamma=2.2;
            [I, R] = jiep(Im);
            hsv = rgb2hsv(Im);
            I_gamma = I.^(1/gamma);
            S_gamma = R .* I_gamma;
            hsv(:,:,3) = S_gamma;
            eIm = hsv2rgb(hsv);
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(eIm);
            LOEs(i) = LOE(eIm, Im);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        elseif strcmp(method, 'WVM_CVPR2016') == 1
            Im=double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            if size(Im,3)>1
                HSV = rgb2hsv(Im);   % RGB space to HSV  space
                S = HSV(:,:,3);       % V layer
            else
                S = Im;              % gray image
            end
            c_1 = 0.01; c_2 = 0.1; lambda = 1;     % set parameters
            epsilon_stop = 1e-3;  % stopping criteria
            [ R, L, epsilon_R, epsilon_L ] = WVM_CVPR2016( S, c_1, c_2, lambda, epsilon_stop );
            %%% Gamma correction
            gamma = 2.2;
            L_gamma = 255*((L/255).^(1/gamma));
            enhanced_V = R .* L_gamma;
            HSV(:,:,3) = enhanced_V;
            eIm = hsv2rgb(HSV);
            imwrite(eIm/255, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(eIm/255);
            LOEs(i) = LOE(eIm/255, Im/255);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        elseif strcmp(method, 'HE') == 1
            Im=im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            [eIm, ~] = histeq(Im);
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(eIm);
            LOEs(i) = LOE(eIm, Im);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        end
        % imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
    end
    matname = [write_mat_dir method '.mat'];
    mNIQEs = mean(NIQEs);
    mLOEs = mean(LOEs);
    fprintf('Mean NIQE = %2.4f, Mean LOE = %2.4f\n', mNIQEs, mLOEs);
    save(matname, 'NIQEs', 'mNIQEs', 'LOEs', 'mLOEs');
end
% subplot(1,2,1); imshow(im); title('Input');
% subplot(1,2,2); imshow(enhance); title('Illumination Adjustment');