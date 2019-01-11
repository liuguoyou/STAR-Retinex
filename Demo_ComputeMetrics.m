clc;
clear;
%%% test dataset
datasets = {'LowLight', 'NASA', 'LDR', 'NPE', 'VV'};
for d = 1:length(datasets)
    Testset = datasets{d}; % select test dataset
    Test_dir  = fullfile('/home/csjunxu/Paper/Enhancement/Dataset', ['Images_' Testset]);
    %%% read images
    ext         =  {'*.jpg','*.jpeg','*.JPG','*.png','*.bmp'};
    im_dir   =  [];
    for i = 1 : length(ext)
        im_dir = cat(1,im_dir, dir(fullfile(Test_dir,ext{i})));
    end
    im_num = length(im_dir);
    
    %%% metrics
    addpath(genpath('metrics'));
    % addpath('/home/csjunxu/Paper/Enhancement/Metrics/vifvec_release');
    metrics = {'LOE', 'NIQE', 'VLD', 'VIF', 'ARISMC', 'AB', 'DE', 'EME', 'PixDist'};
    % LOE: 0~1;
    % NIQE, VLD, VIF, ARISMC: uint8;
    % ARISMC(1) Luminance component only
    % ARISMC(2) Luminance and chromatic components
    
    %%% methods
    addpath('methods');
    methods = {'None', 'JieP_ICCV2017', 'WVM_CVPR2016', 'MF_SP2016', 'SRIE_TIP2015', ...
        'NPE_TIP2013', 'BPDHE_TCE2010', 'MSRCR', 'SSR_TIP1997', 'HE', 'Li_TIP2018'};
    % Li_TIP2018 will run out of memory on 13.bmp
    
    %%% save
    % write_mat_dir  = ['/home/csjunxu/Github/data/Ultrasound/'];
    % write_mat_dir  = ['/home/csjunxu/Github/Segmentation-master/'];
    write_mat_dir = ['/home/csjunxu/Paper/Enhancement/Results_' Testset '/'];
    % write_mat_dir = '/home/csjunxu/Paper/Enhancement/Results_NASA/';
    
    %%% begin
    for m = 1:length(methods)
        method = methods{m};
        if strcmp(method, 'None') == 1
            Enhance_dir = Test_dir;
        else
            Enhance_dir = [write_mat_dir method '/'];
        end
        if ~isdir(Enhance_dir)
            fprintf('No %s Results on %s!\n', method, Testset);
        end
        NIQEs = zeros(im_num,1);
        LOEs = zeros(im_num,1);
        VLDs = zeros(im_num,1);
        VIFs = zeros(im_num,1);
        ARISMCs = zeros(im_num,2);
        ABs = zeros(im_num,1);
        DEs = zeros(im_num,1);
        EMEs = zeros(im_num,1);
        PixDs = zeros(im_num,1);
        eim_dir   =  [];
        for i = 1 : length(ext)
            eim_dir = cat(1,eim_dir, dir(fullfile(Enhance_dir,ext{i})));
        end
        for i = 1:im_num
            Im = imread(fullfile(Test_dir, im_dir(i).name));
            eIm = imread(fullfile(Enhance_dir, eim_dir(i).name));
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
        end
        matname = [write_mat_dir method '.mat'];
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