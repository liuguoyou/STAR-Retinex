clc;clear;
%%% choose test dataset
datasets = {'LowLight', 'NPE', 'VV', 'NASA', 'LDR'};
%%% metrics
addpath(genpath('metrics'));
% addpath('/home/csjunxu/Paper/Enhancement/Metrics/vifvec_release');
metrics = {'NIQE', 'VIF'};
% NIQE and VIF: input image is in uint8;

method = 'STAR';
for d = 1:length(datasets)
    Testset = datasets{d}; % select test dataset
    Test_dir = fullfile('/home/csjunxu/Paper/Enhancement/Dataset', ['Images_' Testset]);
    %%% read images
    ext      =  {'*.jpg','*.jpeg','*.JPG','*.png','*.bmp'};
    im_dir   =  [];
    for i = 1 : length(ext)
        im_dir = cat(1,im_dir, dir(fullfile(Test_dir,ext{i})));
    end
    im_num = length(im_dir);
    write_mat_dir = ['/home/csjunxu/Paper/Enhancement/Results_' Testset '/'];
    write_img_dir = [write_mat_dir method '/'];
    if ~isdir(write_img_dir)
        mkdir(write_img_dir);
    end 
    gamma=2.2;
    alpha = 0.0005;
    beta = 0.0005;
    for pI = [1.5]  
        for pR = [0.5]
            NIQEs = zeros(im_num,1);
            VIFs = zeros(im_num,1);
            for i = 11%1:im_num 
                name = regexp(im_dir(i).name, '\.', 'split');
                Im=im2double( imread(fullfile(Test_dir, im_dir(i).name)) );
                [I, R] = STAR(Im, alpha, beta, pI, pR);
                hsv = rgb2hsv(Im);
                I_gamma = I.^(1/gamma);
                S_gamma = R .* I_gamma;
                hsv(:,:,3) = S_gamma;
                eIm = hsv2rgb(hsv); 
                % convert Im and eIm to uint8
                Im = uint8(Im*255);
                eIm = uint8(eIm*255);
                % metrics 
                NIQEs(i) = niqe(eIm);
                VIFs(i) = VIF(Im,eIm);
                fprintf('%s : NIQE = %2.2f, VIF = %2.2f\n', im_dir(i).name, NIQEs(i), VIFs(i));
                imwrite(eIm, [method '_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
                  num2str(alpha) '_beta=' num2str(beta) '_' name{1} '.png'])
            end
            matname = [write_mat_dir '/Our_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
                num2str(alpha) '_beta=' num2str(beta) '.mat'];
            mNIQEs = mean(NIQEs);
            mVIFs = mean(VIFs);
            fprintf('mNIQE = %2.4f, mVIF = %2.4f\n', mNIQEs, mVIFs);
            save(matname, 'NIQEs', 'mNIQEs', 'VIFs', 'mVIFs');
        end
    end
end