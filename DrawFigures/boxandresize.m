function [ outputimage ] = boxandresize( image, top, left, width, scale, model )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% top ������
% left ���������
% width ����
% scale �Ŵ���
% model ģʽ
outputimage = DrawBox( image, top, left, width, [0,255,0]);

smallimage = outputimage(top:top+width,left:left+width,:);
bigimage = imresize(smallimage,scale, 'bic');
[m ,n, c] = size(bigimage);
if(model==1)    %���½�
    outputimage(end-m+1:end,end-n+1:end,:) = bigimage;
elseif(model==2)   %���½�
    outputimage(end-m+1:end,1:n,:) = bigimage;
else
        outputimage(1:m,end-n+1:end,:) = bigimage;
end

end

