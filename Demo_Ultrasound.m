clc;
clear;
Original_image_dir  =    'dataset/'; % put clean images in this folder
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.jpg');
im_dir  = dir(fpath);
im_num = length(im_dir);

write_dir  = ['/home/csjunxu/Github/data/Ultrasound/'];
if ~isdir(write_dir)
    mkdir(write_dir);
end

gamma=2.2;
alpha = 0.001;
beta = 0.0001;
for i = 1:im_num
    for pI = [.5:.1:3]
        for pR = [.5:.1:3]
            im=im2double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            [I, R] = enhancer(im, alpha, beta, pI, pR);
            hsv = rgb2hsv(im);
            I_gamma = I.^(1/gamma);
            S_gamma = R .* I_gamma;
            hsv(:,:,3) = S_gamma;
            enhance = hsv2rgb(hsv);
            imwrite(enhance, [write_dir im_dir(i).name '_pI=' num2str(pI) '_pR=' num2str(pR) '_alpha=' ...
                num2str(alpha) '_beta=' num2str(beta) '.png']);
        end
    end
end
% subplot(1,2,1); imshow(im); title('Input');
% subplot(1,2,2); imshow(enhance); title('Illumination Adjustment');