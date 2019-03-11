function [ O1, O2] = STAR( src, alpha, beta, pI, pR, vareps, r, K, debug)
if (~exist('alpha','var'))	% alpha -- parameter for structure
    alpha = 0.0005;
end
if (~exist('beta','var'))	% beta -- parameter for texture
    beta = 0.0005;
end
if (~exist('pI','var'))	% pI -- parameter for structure INTENSITY
    pI = 1;
end
if (~exist('pR','var'))	% pR -- parameter for texture INTENSITY
    pR = 1;
end
if (~exist('vareps','var')) % vareps -- stopping parameter
    vareps = 0.01;
end
if (~exist('r','var'))      % r -- the radius of Omega in Eq.(3)
    r = 1;
end
if (~exist('K','var'))      % K -- maximum iterations
    K = 20;
end
if (~exist('debug','var'))  % debug -- set debug/release
    debug = true;
end
eps=0.0001;
if size(src,3)==1
    O = src;
else
    hsv = rgb2hsv(src);
    O = hsv(:,:,3);
end

O1=O.^0.5;   % initialize I_0
O2=O1;        % initialize R_0
if debug == true
    fprintf('-- Stop iteration until eplison < %02f or K > %d\n', vareps, K);
end
for iter = 1:K
    preO1=O1;
    preO2=O2;
    %% algorithm for I
    maxE=max(pI,pR);
    O1=O./O2;
    O1x = diff(O1,1,2); O1x = padarray(O1x, [0 1], 'post');
    O1y = diff(O1,1,1); O1y = padarray(O1y, [1 0], 'post');
    avgO1x=convBox( single(O1x), r);
    avgO1y=convBox( single(O1y), r);
    uO1x = max(abs(avgO1x).^maxE,eps).^(-1);  % structure map avgIx.^pI > avgIx.*Ix > Ix.^2
    uO1y = max(abs(avgO1y).^maxE,eps).^(-1);  % structure map
    uO1x(:,end) = 0;
    uO1y(end,:) = 0;
    
    O1 = solveLinearSystem(O, O2, uO1x, uO1y, alpha);  % Eq.(12)
    epO1 = norm(O1-preO1, 'fro')/norm(preO1, 'fro');   % iterative error of I
    
    %% algorithm for R
    minE=min(pI,pR);
    O2=O./O1;
    O2x = diff(O2,1,2); O2x = padarray(O2x, [0 1], 'post');
    O2y = diff(O2,1,1); O2y = padarray(O2y, [1 0], 'post');
    avgO2x=convBox( single(O2x), r);
    avgO2y=convBox( single(O2y), r);
    vO2x = max(abs(avgO2x).^minE,eps).^(-1);  % texture map
    vO2y = max(abs(avgO2y).^minE,eps).^(-1);  % texture map
    vO2x(:,end) = 0;
    vO2y(end,:) = 0;
    
    O2 = solveLinearSystem(O, O1, vO2x, vO2y, beta);            	% Eq.(13)
    epO2 = norm(O2-preO2, 'fro')/norm(preO2, 'fro');   % iterative error of R
    %% iteration until convergence
    if debug == true
        fprintf('Iter #%d : eplisonI = %f; eplisonR = %f\n', iter, epO1, epO2);
    end
    if(epO1<vareps||epO2<vareps)
        break;
    end
end
O1(O1<0)=0;
O2(O2<0)=0;
end


