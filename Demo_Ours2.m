clc;clear;
%%% choose test dataset
datasets = {'NPE', 'LowLight', 'NASA', 'LDR', 'VV'};
%%% metrics
addpath(genpath('metrics'));
% addpath('/home/csjunxu/Paper/Enhancement/Metrics/vifvec_release');
metrics = {'LOE', 'NIQE', 'VLD', 'VIF', 'ARISMC', 'AB', 'DE', 'EME', 'PixDist'};
% LOE: 0~1;
% NIQE, VLD, VIF, ARISMC: uint8;
% ARISMC(1) Luminance component only
% ARISMC(2) Luminance and chromatic components

method = 'Our';
for d = 2%1:length(datasets)
    Testset = datasets{d}; % select test dataset
    Test_dir  = fullfile('/home/csjunxu/Paper/Enhancement/Dataset', ['Images_' Testset]);
    %%% read images
    ext         =  {'*.jpg','*.jpeg','*.JPG','*.png','*.bmp'};
    im_dir   =  [];
    for i = 1 : length(ext)
        im_dir = cat(1,im_dir, dir(fullfile(Test_dir,ext{i})));
    end
    im_num = length(im_dir);
    % write_mat_dir  = ['/home/csjunxu/Github/data/Ultrasound/'];
    % write_mat_dir  = ['/home/csjunxu/Github/Segmentation-master/'];
    write_mat_dir = ['/home/csjunxu/Paper/Enhancement/Results_' Testset '/'];
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
            VLDs = zeros(im_num,1);
            VIFs = zeros(im_num,1);
            ARISMCs = zeros(im_num,2);
            ABs = zeros(im_num,1);
            DEs = zeros(im_num,1);
            EMEs = zeros(im_num,1);
            PixDs = zeros(im_num,1);
            for i = 1:im_num
                name = regexp(im_dir(i).name, '\.', 'split');
                Im=im2double( imread(fullfile(Test_dir, im_dir(i).name)) );
                [I, R] = enhancer(Im, alpha, beta, pI, pR);
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
                LOEs(i) = LOE(im2double(eIm), im2double(Im));
                VLDs(i) = VLD(eIm, Im);
                VIFs(i) = VIF(Im,eIm);
                ARISMCs(i,:) = ARISMC(eIm);
                ABs(i) = AB(eIm);
                DEs(i) = DE(eIm);
                EMEs(i) = EME(double(eIm));
                PixDs(i) = PixDist(eIm);
                fprintf('%s : NIQE = %2.2f, LOE = %2.2f, VLD = %2.2f, VIF = %2.2f, AB = %2.2f, DE = %2.2f\n', ...
                    im_dir(i).name, NIQEs(i), LOEs(i), VLDs(i), VIFs(i), ABs(i), ...
                    DEs(i));
                % imwrite(enhance, [write_img_dir method '_' name{1} '_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
                %    num2str(alpha) '_beta=' num2str(beta) '_' name{1} '.jpg'])
            end
            matname = [write_mat_dir '/Our_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
                num2str(alpha) '_beta=' num2str(beta) '.mat'];
            mNIQEs = mean(NIQEs);
            mLOEs = mean(LOEs);
            mVLDs = mean(VLDs);
            mVIFs = mean(VIFs);
            mARISMCs = mean(ARISMCs);
            mABs = mean(ABs);
            mDEs = mean(DEs);
            mEMEs = mean(EMEs);
            mPixDs = mean(PixDs);
            fprintf('mNIQE = %2.4f, mLOE = %2.4f, mVLD = %2.4f, mVIF = %2.4f\n', ...
                mNIQEs, mLOEs, mVLDs, mVIFs, mABs);
            save(matname, 'NIQEs', 'mNIQEs', 'LOEs', 'mLOEs', 'VLDs', 'mVLDs', ...
                'VIFs', 'mVIFs', 'ARISMCs', 'mARISMCs', 'ABs', 'mABs', ...
                'DEs', 'mDEs', 'EMEs', 'mEMEs', 'PixDs', 'mPixDs');
        end
    end
end