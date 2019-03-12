clc;
clear;
close all;

I = im2double(imread('correction/img/books.png'));
for scale=1.4
    for alpha=0.001
        for beta=0.0001
            for pI=2
                for pR=1.2:0.1:1.5
                    for ch=1:size(I,3)
                        [L(:,:,ch),R(:,:,ch)]=STAR(I(:,:,ch),alpha,beta,pI,pR);
                    end
                    r=mean(reshape(R,size(R,1)*size(R,2),size(R,3)));
                    for i=1:3
                        temp=L(:,:,i);
                        meanLi=mean(temp(:));
                        L2(:,:,i)=L(:,:,i)./meanLi.*(scale*mean(L(:))*r(i));
                    end
                    I_c=L2.*R;
                    %imshow(R);
                    %imshow(L2);
                    %imshow(I_c);
                    imwrite(I_c,['/home/csjunxu/Paper/Enhancement/Results_Color/books_star_' num2str(scale) '_' num2str(alpha) '_' num2str(beta) '_' num2str(pI) '_' num2str(pR) '.png']);
                end
            end
        end
    end
end
