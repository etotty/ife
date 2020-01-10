

function [betaIFEH, betaiIFE, varcovIFEi, seIFEi, sigma2, SSR, nnn, r1, r2, betaiIFE0, betaiOLS, F1, L1, F2, L2, VNT1, VNT2] = IFEifunctionNEW(X1,Y1,order,S,rmax,rfix, betaIFEp, seIFEp, mgfe); 
[T,N,p]=size(X1);



% remove state-specific trends and person fixed effects
[Xdet,Ydet] = RemoveTimeTrends(X1,Y1,order,S);

if mgfe==1
[X,Y]=TimeDemean(Xdet,Ydet,1);
elseif mgfe==0
Y=Ydet;
X=Xdet;
end

% calculate XXinv outside of loop to save time
[XXinv]= Mul_XXinv(X); 




% Estimation
% betaIFEH=zeros(p,2);
betaiIFE0=zeros(p,2*N);
betaiIFE=zeros(p,2*N);    % to contain the interative effect estimators, with different staring methods (2 methods)
nnn=zeros(1,2);   % contain the number of iterations for the 2 methods to achieve convergence
sigma2=zeros(1,2); % the estimated residual variance for each starting method, also the optimal value of the objective function
seIFEi=zeros(p,2);     % contains the IFE standard errors for the 2 different starting methods
SSR=zeros(1,2);
betaIFEH=zeros(p,2);


% Starting method 1:
% Initiate IFE by calculating standard panel beta estimator (lambda=0) as a starting value
[betaiOLS]=Mul_panelbetai(X,Y); 

U=Y;
for i=1:N
for k=1:p;
    U(:,i)=U(:,i)-X(:,i,k)* betaiOLS(k,i);
end
end

% Get the first estimate of r and F,lambda (because iteration below starts by estimating beta
if rfix==0
    [r1,~,~,~]=nbplog(U,rmax,1,0);
elseif rfix==1
    r1=rmax;
end
[F1,L1,VNT]=panelFactorNew(U,r1);

% Perform IFE iteration
[betaiIFE0(:,1:N), F1,L1, VNT1, e1, nnn(1,1),r1]=Mul_betaiIterNewET(X,XXinv, Y, F1,L1, r1,rmax,rfix, 0.00000001, betaIFEp, seIFEp);
sigma2(1,1)=trace(e1*e1')/(N*T-r1*(N+T)+r1^2-2);
SSR(1,1)=trace(e1*e1');
betaIFEH(:,1)=mean(betaiIFE0(:,1:N),2);
% seIFE(:,1)=seBAIET(X,F1,L1,e1);
squaredIFEi=zeros(p,p,N);
for i=1:N
    squaredIFEi(:,:,i)=(betaiIFE0(:,i)-betaIFEH(:,1))*(betaiIFE0(:,i)-betaIFEH(:,1))';
end;
varcovIFEi=(1/(N*(N-1)))*sum(squaredIFEi,3);
seIFEi(:,1)=sqrt(diag(varcovIFEi));

% correct beta for serial corr, cross-sec corr, and heteroskedasticity
% betaiIFE(:,1:N)=biasiBAIET(X,F1,L1,e1,betaiIFE0(:,1:N));



%Starting method 2:
% Initiate IFE by calculating the factors and factor loading (beta=0)
if rfix==0
    [r2,~,~,~]=nbplog(Y,rmax,1,0);
elseif rfix==1
    r2=rmax;
end
[F12,L12,VNT2]=panelFactorNew(Y,r2);

%Perform IFE iteration
[betaiIFE0(:,N+1:2*N), F2,L2, VNT2, e2, nnn(1,2),r2]=Mul_betaiIterNewET(X,XXinv, Y, F12,L12, r2,rmax,rfix, 0.00000001, betaIFEp, seIFEp);
sigma2(1,2)=trace(e2*e2')/(N*T-r2*(N+T)+r2^2-2);
SSR(1,2)=trace(e2*e2');
% seIFE(:,2)=seBAIET(X,F2,L2,e2);
betaIFEH(:,2)=mean(betaiIFE0(:,N+1:2*N),2);
squaredIFEi=zeros(p,p,N);
for i=1:N
    squaredIFEi(:,:,i)=(betaiIFE0(:,N+i)-betaIFEH(:,2))*(betaiIFE0(:,N+i)-betaIFEH(:,2))';
end;
varcovIFEi=(1/(N*(N-1)))*sum(squaredIFEi,3);
seIFEi(:,2)=sqrt(diag(varcovIFEi));

% correct beta for serial corr, cross-sec corr, and heteroskedasticity
% betaiIFE(:,N+1:2*N)=biasiBAIET(X,F2,L2,e2,betaiIFE0(:,N+1:2*N));




% bootstrap for significance
% if sigma2(1,1)<=sigma2(1,2)
%     for 
% [tstat,cvBOOT]=WildClusterBootstrapT(Y,X,betaIFE(:,1),seIFE(:,1),F1,L1,XXinv);
% else %send method 2 results in for teen, because same ssr and I preer method 2 results
% [tstat,cvBOOT]=WildClusterBootstrapT(Y,X,betaIFE(:,2),seIFE(:,2),F2,L2,XXinv);    


betaiIFE;
seIFEi;
sigma2;
% betaIFEH(:,1)=mean(betaiIFE(:,1:N),2);
% betaIFEH(:,2)=mean(betaiIFE(:,N+1:2*N),2);



