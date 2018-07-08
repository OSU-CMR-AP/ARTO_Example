function [phi_Cor,cMap] = wls(phi,sigma,M,opt)

% Build A matrix
[A,Af] = getMatA(M,opt.pOrd);

% Weight data
tmp = sigma(M);

A_w = bsxfun(@times, A, 1./tmp);
phi_w = phi(M).*(1./tmp);
    
AHAm = A_w'*A_w;

AHy = A_w'*phi_w;

cLS   = AHAm \ AHy;                 % Coefficients
cMap = reshape(Af*cLS,size(M));     % fitted maps 

phi_Cor = phi - cMap;
    
end

