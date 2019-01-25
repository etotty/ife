% Iterates between estimating factor structure given regressions coefficients, and vice-versa, until convergence. This code was kindly provided to me by Bia via email correspondence.

% OUTPUT:
% betanew: computed beta under iteration with error precision=tolerate
% factor: estimated factor
% lambda: estimated loadings
% V: the eigenvalues matrix
% e: estimated residuals
% niter: number of interations to achieve convergence 

function [betanew, factor, lambda, V, e, niter,r]=Mul_betaIterNew(X, xxinv, Y, F,L, r,rmax,rfix, tolerate);
   [T,N,p]=size(X);
   changeU2=1;
   sumU2old=.0000000001;

n=0;   % number of iterations needed to stop
while (changeU2  > tolerate & n < 500) 
    n=n+1;
    [beta]=Mul_panelbetaNew(X,xxinv, Y, F,L);
    U=Y;
    for k=1:p;
       U=U-X(:,:,k)* beta(k,1);
    end
    U2=U.^2;
    sumU2=sum(U2(:));
    changeU2=(sumU2-sumU2old)/sumU2old;
    sumU2old=sumU2;

    if rfix==0
      [r,~,~,~]=nbplog(U,rmax,1,0);
    end
    [F,L, VNT]=panelFactorNew(U,r);
  
    
end
betanew=beta;
niter=n;
factor=F;
lambda=L;
V=VNT;
e=U-F*L'; 