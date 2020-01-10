% OUTPUT:
% betanew: computed beta under iteration with error precision=tolerate
% factor: estimated factor
% lambda: estimated loadings
% V: the eigenvalues matrix
% e: estimated residuals
% niter: number of interations to achieve convergence 

function [betainew, factor, lambda, V, e, niter,r]=Mul_betaiIterNewET(X, xxinv, Y, F,L, r,rmax,rfix, tolerate, betaIFEp, seIFEp);
   [T,N,p]=size(X);
   changeU2=1;
   sumU2old=.0000000001;

%     tolerate=seIFEp*2.575;
%     maxgap=seIFEp*10;

% maxsqrtsumBetadiffsquared=1;

n=0;   % number of iterations needed to stop
while (changeU2  > tolerate & n < 300) 
    n=n+1;
    [betai]=Mul_panelbetaiNew(X,xxinv, Y, F,L);
    U=Y;
    for i=1:N
    for k=1:p;
       U(:,i)=U(:,i)-X(:,i,k)* betai(k,i);
    end
    end
    
    U2=U.^2;
    sumU2=sum(U2(:));
    changeU2=(sumU2-sumU2old)/sumU2old;
    sumU2old=sumU2;

%     gap=abs(betai-betaIFEp);
%     maxgap=max(gap);

%     Betadiff=zeros(p,N);
%     for i=1:N
%     for k=1:p
%         Betadiff(k,i)=betai(k,i)-betaIFEp(k,1);
%     end
%     end
%     Betadiffsquared=Betadiff.^2;
%     sumBetadiffsquared=sum(Betadiffsquared,2);
%     sqrtsumBetadiffsquared=sqrt(sumBetadiffsquared);
%     maxsqrtsumBetadiffsquared=max(sqrtsumBetadiffsquared);
        

    if rfix==0
      [r,~,~,~]=nbplog(U,rmax,1,0);
    end
    [F,L, VNT]=panelFactorNew(U,r);
  
    
end
betainew=betai;
niter=n;
factor=F;
lambda=L;
V=VNT;
e=U-F*L'; 