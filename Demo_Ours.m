clc;
clear;
Original_image_dir  =    'dataset/'; % put clean images in this folder
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.bmp');
im_dir  = dir(fpath);
im_num = length(im_dir);

% write_dir  = ['/home/csjunxu/Github/data/Ultrasound/'];
% write_dir  = ['/home/csjunxu/Github/Segmentation-master/'];
write_dir  = ['./NIQEResults/'];
if ~isdir(write_dir)
    mkdir(write_dir);
end

gamma=2.2;
alpha = 0.001;
beta = 0.0001;
for pI = [.1:.1:2]
    for pR = [.1:.1:2]
        score = zeros(im_num,1);
        for i = 1:im_num
            Im=im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            [I, R] = enhancer(Im, alpha, beta, pI, pR);
            hsv = rgb2hsv(Im);
            I_gamma = I.^(1/gamma);
            S_gamma = R .* I_gamma;
            hsv(:,:,3) = S_gamma;
            eIm = hsv2rgb(hsv);
            score(i) = niqe(eIm);
            fprintf('%s : NIQE = %2.4f\n', im_dir(i).name, score(i));
            % imwrite(enhance, [write_dir im_dir(i).name '_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
            %    num2str(alpha) '_beta=' num2str(beta) '.jpg'])
        end
        matname = ['NIQEResults/Our_aIpI=' num2str(pI) '_RpR=' num2str(pR) '_alpha=' ...
            num2str(alpha) '_beta=' num2str(beta) '.mat'];
        mscore = mean(score);
        if mscore < 2.80
            fprintf('Mean NIQE = %2.4f\n', mscore);
            save(matname, 'score', 'mscore');
        end
    end
end
% subplot(1,2,1); imshow(im); title('Input');
% subplot(1,2,2); imshow(enhance); title('Illumination Adjustment');