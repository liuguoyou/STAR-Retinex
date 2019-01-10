% Single Scale Retinex, TIP 1997


%clear,clc;
%Img = imread('/home/csjunxu/Github/JieP/dataset/campus.bmp');
%eImg = SSR_TIP1997(Img, 10000);

%subplot(1,2,1); imshow(Img); title('Input');
%subplot(1,2,2); imshow(eImg); title('SSR Enhanced Adjustment');




function eImg = SSR_TIP1997(Img, c)
[h, w, ch] = size(Img);
sum = 0;
for i=1:h
    for j=1:w
        sum=sum+exp(-(i^2+j^2)/c^2);
    end
end
k=1/sum;
for t=1:3
    wei=double(Img(:,:,t));
    for i=1:h
        for j=1:w
            f(i,j)=k*exp(-(i^2+j^2)/c^2);
        end
    end
    q=double(conv2(wei,f,'same'));
    eImg(:,:,t)=log(wei)-log(q);
end
end
    