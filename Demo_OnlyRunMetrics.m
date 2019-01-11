clc;
clear;
% test images
Original_image_dir  =    '/home/csjunxu/Paper/Enhancement/Dataset/LowLightImages/'; 
% Original_image_dir  =    '/home/csjunxu/Paper/Enhancement/Dataset/NASA/'; 
% Original_image_dir  =    '/home/csjunxu/Paper/Enhancement/Dataset/LDR_TEST_IMAGES_DICM/'; 
fpath = fullfile(Original_image_dir, '*.jpg');
im_dir  = dir(fpath);
im_num = length(im_dir);

% metrics
addpath('metrics');
addpath('/home/csjunxu/Paper/Enhancement/Metrics/vifvec_release');
metrics = {'LOE', 'NIQE', 'VLD'};
% LOE:0~1; NIQE uint8
% methods
addpath('methods');
methods = {'JieP_ICCV2017', 'WVM_CVPR2016', 'MF_SP2016', 'SRIE_TIP2015', ...
    'NPE_TIP2013', 'BPDHE_TCE2010', 'MSRCR', 'SSR_TIP1997', 'HE', 'Li_TIP2018'};
% Li_TIP2018 will run out of memory on 13.bmp

for m = 6%1:length(methods)
    method = methods{m};
    % write_mat_dir  = ['/home/csjunxu/Github/data/Ultrasound/'];
    % write_mat_dir  = ['/home/csjunxu/Github/Segmentation-master/'];
    % write_mat_dir = '/home/csjunxu/Paper/Enhancement/Results_LowLight/';
    % write_mat_dir = '/home/csjunxu/Paper/Enhancement/Results_NASA/';
    write_mat_dir = '/home/csjunxu/Paper/Enhancement/Results_LDR/';
    write_img_dir = [write_mat_dir method '/'];
    if ~isdir(write_img_dir)
        mkdir(write_img_dir);
    end
    NIQEs = zeros(im_num,1);
    LOEs = zeros(im_num,1);
    VLDs = zeros(im_num,1);
    for i = 1:im_num
        name = regexp(im_dir(i).name, '\.', 'split');
        if strcmp(method, 'BPDHE_TCE2010') == 1
            Enhance_image_dir = ['/home/csjunxu/Paper/Enhancement/Results_LowLight/' method]; 
            efpath = fullfile(Enhance_image_dir, '*.jpg');
            eim_dir  = dir(efpath);
            Im = imread(fullfile(Original_image_dir, im_dir(i).name));
            eIm = imread(fullfile(Enhance_image_dir, eim_dir(i).name));
            NIQEs(i) = niqe(eIm);
            LOEs(i) = LOE(im2double(eIm), im2double(Im));
            VLDs(i) = VLD(eIm, Im);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, ...
                NIQEs(i), LOEs(i));
        elseif strcmp(method, 'MSRCR') == 1
            addpath('methods/multiscaleRetinex/');
            Im = imread(fullfile(Original_image_dir, im_dir(i).name));
            eIm = multiscaleRetinex(Im, 'MSRCR');
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(uint8(eIm));
            LOEs(i) = LOE(eIm/255, Im/255);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, ...
                NIQEs(i), LOEs(i));
        elseif strcmp(method, 'SSR_TIP1997') == 1
            Im = imread(fullfile(Original_image_dir, im_dir(i).name));
            eIm = SSR_TIP1997(Im, 10000);
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(uint8(eIm));
            LOEs(i) = LOE(eIm/255, Im/255);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, ...
                NIQEs(i), LOEs(i));
        elseif strcmp(method, 'NPE_TIP2013') == 1
            addpath('methods/NPE_TIP2013/');
            Im = imread(fullfile(Original_image_dir, im_dir(i).name));
            eIm=NPEA(fullfile(Original_image_dir, im_dir(i).name));
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(eIm);
            LOEs(i) = LOE(im2double(eIm), im2double(Im));
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, ...
                NIQEs(i), LOEs(i));
        elseif strcmp(method, 'SRIE_TIP2015') == 1
            Im = imread(fullfile(Original_image_dir, im_dir(i).name));
            alpha = 1000; beta= 0.01; gamma = 0.1; lambda = 10; % set parameters
            error_R = 10; error_I = 10; % initial stopping criteria error_R and error_I
            stop = 0.1;  % stopping criteria
            HSV = rgb2hsv( double(Im) );   % RGB space to HSV  space
            S = HSV(:,:,3);       % V layer
            [ R, I, error_R, error_I ] = SRIE_TIP2015( S, alpha, beta, gamma, lambda, ...
                error_R, error_I, stop);
            % Gamma correction
            gamma1 = 2.2;
            I_gamma = 255 * ( (I/255).^(1/gamma1) );
            enhanced_V = R .* I_gamma;
            HSV(:,:,3) = enhanced_V;
            eIm = hsv2rgb(HSV);  %  HSV space to RGB space
            eIm = cast(eIm, 'uint8');
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(eIm);
            LOEs(i) = LOE(im2double(eIm), im2double(Im));
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
            NIQEs(i) = niqe(uint8(eIm*255));
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
            NIQEs(i) = niqe(uint8(eIm));
            LOEs(i) = LOE(eIm/255, Im/255);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        elseif strcmp(method, 'MF_SP2016') == 1
            Im=double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            eIm = MF_SP2016(Im);
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(uint8(eIm));
            LOEs(i) = LOE(eIm/255, Im/255);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        elseif strcmp(method, 'HE') == 1
            Im=im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            [eIm, ~] = histeq(Im);
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(uint8(eIm*255));
            LOEs(i) = LOE(eIm, Im);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        elseif strcmp(method, 'Li_TIP2018') == 1
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
            imwrite(eIm, [write_img_dir method '_' name{1} '.jpg']);
            NIQEs(i) = niqe(uint8(eIm));
            LOEs(i) = LOE(eIm/255, Im/255);
            fprintf('%s : NIQE = %2.4f, LOE = %2.4f\n', im_dir(i).name, NIQEs(i), LOEs(i));
        end
    end
    matname = [write_mat_dir method '.mat'];
    mNIQEs = mean(NIQEs);
    mLOEs = mean(LOEs);
    fprintf('Mean NIQE = %2.4f, Mean LOE = %2.4f\n', mNIQEs, mLOEs);
    save(matname, 'NIQEs', 'mNIQEs', 'LOEs', 'mLOEs');
end
% subplot(1,2,1); imshow(im); title('Input');
% subplot(1,2,2); imshow(enhance); title('Illumination Adjustment');