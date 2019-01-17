clc;clear;
%%% choose test dataset
datasets = {'BSDS_LL'};
ext         =  {'*.jpg','*.jpeg','*.JPG','*.png','*.bmp'};
nSig = 5; % noise level of additive white Gaussian noise
for d = 1:length(datasets)
    Testset = datasets{d}; % select test dataset
    Test_dir  = fullfile('/home/csjunxu/Dataset', Testset, 'LowLightImages');
    Original_dir  = fullfile('/home/csjunxu/Dataset', Testset, 'OriginalImages');
    %%% read images
    oIm_dir   =  [];
    for j = 1 : length(ext)
        oIm_dir = cat(1,oIm_dir, dir(fullfile(Original_dir,ext{j})));
    end
    im_dir   =  [];
    for j = 1 : length(ext)
        im_dir = cat(1,im_dir, dir(fullfile(Test_dir,ext{j})));
    end
    im_num = length(im_dir);
    %%% methods
    addpath(genpath('../'));
    methods = {'LIME_TIP2017', 'JieP_ICCV2017', 'WVM_CVPR2016', 'MF_SP2016', ...
        'NPE_TIP2013', 'SRIE_TIP2015', 'LDR_TIP2013', 'CVC_TIP2011', ...
        'WAHE_TIP2009', 'BPDHE_TCE2010', 'MSRCR', 'SSR_TIP1997', 'HE', ...
        'Dong_ICME2011', 'BIMEF_2019', 'None', 'Li_TIP2018'};
    % 'Li_TIP2018': run out of memory or SVD include NaN or Inf

    %%% begin comparisons
    write_mat_dir = ['/home/csjunxu/Dataset/' Testset '/'];
    % write_mat_dir = '/home/csjunxu/Paper/Enhancement/Results_NASA/';
    for m = 1:length(methods)
        method = methods{m};
        write_img_dir = [write_mat_dir method '/'];
        if ~isdir(write_img_dir)
            mkdir(write_img_dir);
        end
        % record all the results in each iteration
        PSNRs = zeros(im_num, 1, 'double');
        SSIMs = zeros(im_num, 1, 'double');
        for i = 1:im_num
            name = regexp(im_dir(i).name, '\.', 'split');
            oIm =  im2double(imread(fullfile(Original_dir, oIm_dir(i).name))); % original image
            lIm =  im2double( imread(fullfile(Test_dir, im_dir(i).name)) ); % low-light image
            randn('seed',0);
            lnIm = lIm + nSig/255*randn(size(lIm)); % low-light noisy image
            lnIm(lnIm<0)=0; % avoid negative values, may generate bugs
            PSNR = csnr( lnIm*255, oIm*255, 0, 0 );
            SSIM = cal_ssim( lnIm*255, oIm*255, 0, 0 );
            fprintf('%s: initial value of PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name,PSNR,SSIM);
            if strcmp(method, 'LIME_TIP2017') == 1
                %addpath(genpath('../methods/BM3D/'));
                Im = lnIm;
                %--------------------------------------------------------------
                post = false; % Denoising?
                para.lambda = .15; % Trade-off coefficient
                % Although this parameter can perform well in a relatively large range,
                % it should be tuned for different solvers and weighting strategies due to
                % their difference in value scale.
                % Typically, lambda for exact solver < for sped-up solver
                % and using Strategy III < II < I
                % ---> lambda = 0.15 is fine for SPED-UP SOLVER + STRATEGY III
                % ......
                para.sigma = 2; % Sigma for Strategy III
                para.gamma = 0.7; %  Gamma Transformation on Illumination Map
                para.solver = 1; % 1: Sped-up Solver; 2: Exact Solver
                para.strategy = 3;% 1: Strategy I; 2: II; 3: III
                %---------------------------------------------------------------
                [eIm, T_ini,T_ref] = LIME_TIP2017(Im,para);
                %% Post Processing
                if post
                    YUV = rgb2ycbcr(eIm);
                    Y = YUV(:,:,1);
                    
                    sigma_BM3D = 10;
                    [~, Y_d] = BM3D(Y,Y,sigma_BM3D,'lc',0);
                    
                    I_d = ycbcr2rgb(cat(3,Y_d,YUV(:,:,2:3)));
                    I_f = (eIm).*repmat(T_ref,[1,1,3])+I_d.*repmat(1-T_ref,[1,1,3]);
                end
                % convert eIm to uint8
                eIm = uint8(eIm*255);
            elseif strcmp(method, 'Dong_ICME2011') == 1
                Im=uint8(lnIm*255);
                eIm = Dong_ICME2011(Im);
                % convert eIm to uint8
                eIm = uint8(eIm*255);
            elseif strcmp(method, 'JieP_ICCV2017') == 1
                Im=lnIm;
                gamma=2.2;
                [I, R] = jiep(Im);
                hsv = rgb2hsv(Im);
                I_gamma = I.^(1/gamma);
                S_gamma = R .* I_gamma;
                hsv(:,:,3) = S_gamma;
                eIm = hsv2rgb(hsv);
                % convert eIm to uint8
                eIm = uint8(eIm*255);
            elseif strcmp(method, 'WVM_CVPR2016') == 1
                Im=lnIm*255;
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
                % convert eIm to uint8
                Im  = uint8(Im);
                eIm = uint8(eIm);
            elseif strcmp(method, 'MF_SP2016') == 1
                Im  = lnIm*255;
                eIm = MF_SP2016(Im);
                % convert eIm to uint8
                Im = uint8(Im);
                eIm = uint8(eIm);
            elseif strcmp(method, 'SRIE_TIP2015') == 1
                Im = uint8(lnIm*255);
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
                % convert eIm to uint8
                eIm = uint8(eIm);
            elseif strcmp(method, 'NPE_TIP2013') == 1
                %addpath('methods/NPE_TIP2013/');
                Im = uint8(lnIm*255);
                eIm=NPEA(fullfile(Test_dir, im_dir(i).name));
            elseif strcmp(method, 'BPDHE_TCE2010') == 1
                Im = uint8(lnIm*255);
                eIm = BPDHE_TCE2010(Im);
            elseif   strcmp(method, 'MSRCR') == 1
                %addpath('methods/multiscaleRetinex/');
                Im = uint8(lnIm*255);
                eIm = multiscaleRetinex(Im, 'MSRCR');
                % convert eIm to uint8
                eIm = uint8(eIm*255);
            elseif strcmp(method, 'SSR_TIP1997') == 1
                Im = uint8(lnIm*255);
                if size(Im,3)>1
                    HSV = rgb2hsv(Im);   % RGB space to HSV  space
                    X = HSV(:,:,3);       % V layer
                else
                    X = Im;              % gray image
                end
                [R,L] = SSR_TIP1997(X, [], 1);
                %%% Gamma correction
                gamma = 2.2;
                L_gamma = ((L/255).^(1/gamma));
                enhanced_V = R .* L_gamma;
                HSV(:,:,3) = enhanced_V;
                eIm = hsv2rgb(HSV);
                % convert eIm to uint8
                eIm = uint8(eIm);
                % eIm = multiscaleRetinex(Im, 'SSR');
                % eIm = SSR_TIP1997(Im, 10000); % will fail on NPE dataset
            elseif strcmp(method, 'HE') == 1
                Im=lnIm;
                [eIm, ~] = histeq(Im);
                % convert eIm to uint8
                Im = uint8(Im*255);
                eIm = uint8(eIm*255);
            elseif strcmp(method, 'BIMEF_2019') == 1
                %addpath('methods/BIMEFutil/');
                Im=uint8(lnIm*255);
                eIm = BIMEF_2019(Im);
                % convert eIm to uint8
                eIm = uint8(eIm);
            elseif strcmp(method, 'LDR_TIP2013') == 1
                %addpath('methods/LDR_TIP2013/');
                Im=uint8(lnIm*255);
                eIm = LDR_TIP2013(Im, 2.5);
            elseif strcmp(method, 'CVC_TIP2011') == 1
                Im=uint8(lnIm*255);
                eIm = CVC_TIP2011(Im);
            elseif strcmp(method, 'WAHE_TIP2009') == 1
                Im=uint8(lnIm*255);
                eIm = WAHE_TIP2009(Im,1.5);
            elseif strcmp(method, 'Li_TIP2018') == 1
                Im = lnIm*255;
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
                [R, L, N] = Li_TIP2018(Im, para); % N is the noise map
                % imshow((N-min(min(N)))./(max(max(N))-min(min(N))))
                eIm = R.*L.^(1/gamma);
                % convert eIm to uint8
                eIm = uint8(eIm);
            end
            %%% image level metrics
            fprintf([Testset ', ' method ', ' name{1} ' is done\n']);
            imwrite(eIm, [write_img_dir method '_n' num2str(nSig) '_' name{1} '.jpg']);
            oIm = uint8(oIm*255);
            PSNRs(i) = csnr( eIm, oIm, 0, 0 );
            SSIMs(i) = cal_ssim( eIm, oIm, 0, 0  );
            fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n',im_dir(i).name,PSNRs(i),SSIMs(i));
        end
        mPSNR=mean(PSNRs);
        mSSIM=mean(SSIMs);
        fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR,mSSIM);
        name = sprintf([write_mat_dir method '_n' num2str(nSig) '.mat']);
        save(name,'nSig','PSNR','SSIM','mPSNR','mSSIM');
    end
end