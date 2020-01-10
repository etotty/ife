%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produces Song's (2013) heterogeneous coefficients version of 
% Bai's (2009) interactive fixed effects (IFE) estimator.
% Bai, Jushan. 2009. "Panel Data Models With Interactive Fixed Effects." Econometrica, 77 (4), 1229-1279.
% Song, M. 2013. "Asymptotic theory for dynamic heterogeneous panels with cross-sectional dependence and its applications." Unpublished manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUTS:
% X1 = (TxNxp) matrix of independent variables
% Y1 = (TxN) matrix dependent variable
% order = the order of a unit-specic trend to be removed from data. order=0 specifies that no trend is removed
% S = matrix that places each cross-section, i, unit into one of M groups
%	(M<N) (i.e., grouping each county into a specific state). This is used for
%	removing group rather than unit trends. Set S=0 if no grouping is
%	necessary
% rmax = maximum number of factors allowed, when using information criterion to select number 
%	of factors (Bai and Ng, 2002)
% rfix = fixes the number of factors at rmax, rather than estimating the number via information 
%	criterion
% betaIFEp = (px1) vector of pooled IFE regression coefficients 
% seIFEp = (px1) vector of standard errors for pooled IFE regression coefficients 
% mgfe = indicator for whether unit-specific regressions should include an intercept (1-yes, 0-no)

% OUTPUTS:
% betaIFEH = (px2) vector of mean-group IFE regression coefficients for each starting method
% betaiIFE = (px2*N) vector of individual IFE regression coefficients for each unit and each starting method
% seIFEi = (px2) vector of standard errors for IFE regression coefficients for each starting method
% sigma2 = (1,2) vector of estimated residual variance for each starting method
% SSR = (1,2) vector of final sum of squared residuals for each starting method
% nnn = (1,2) vector of number of iterations before convergence for each starting method
% r1 = final estimate of number of factors for starting method 1
% r2 = final estimate of number of factors for starting method 2
% F1 = estimated factors for starting method 1
% L1 = estimated loadings for starting method 1
% F2 = estimated factors for starting method 2
% L2 = estiamted loadings for starting method 2
% VNT1 = diagnonal matrix with the r-largest eigenvalues of the second moment matrix of the regression residuals
% VNT2 = diagnonal matrix with the r-largest eigenvalues of the second moment matrix of the regression residuals



function [betaIFEH, betaiIFE, seIFEi, sigma2, SSR, nnn, r1, r2, F1, L1, F2, L2, VNT1, VNT2] = IFEifunction(X1,Y1,order,S,rmax,rfix, betaIFEp, seIFEp, mgfe); 

[T,N,p]=size(X1);


%%%%% 1. Remove state-specific trends and demean
[Xdet,Ydet] = RemoveTimeTrends(X1,Y1,order,S);
if mgfe==1
[X,Y]=TimeDemean(Xdet,Ydet,1);
elseif mgfe==0
Y=Ydet;
X=Xdet;
end


%%%%% 2. Estimation
% calculate XXinv outside of loop to save time
[XXinv]= xxinv(X); 

betaIFEH=zeros(p,2);
betaiIFE=zeros(p,2*N);    % to contain the interative effect estimators, with different staring methods (2 methods)
nnn=zeros(1,2);   % contain the number of iterations for the 2 methods to achieve convergence
sigma2=zeros(1,2); % the estimated residual variance for each starting method, also the optimal value of the objective function
seIFEi=zeros(p,2);     % contains the IFE standard errors for the 2 different starting methods
SSR=zeros(1,2);


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
[betaiIFE(:,1:N), F1,L1, VNT1, e1, nnn(1,1),r1]=Mul_betaiIterNew(X,XXinv, Y, F1,L1, r1,rmax,rfix, 0.00000001, betaIFEp, seIFEp);
sigma2(1,1)=trace(e1*e1')/(N*T-r1*(N+T)+r1^2-2);
SSR(1,1)=trace(e1*e1');
betaIFEH(:,1)=mean(betaiIFE0(:,1:N),2);
squaredIFEi=zeros(p,p,N);
for i=1:N
    squaredIFEi(:,:,i)=(betaiIFE0(:,i)-betaIFEH(:,1))*(betaiIFE0(:,i)-betaIFEH(:,1))';
end;
varcovIFEi=(1/(N*(N-1)))*sum(squaredIFEi,3);
seIFEi(:,1)=sqrt(diag(varcovIFEi));



%Starting method 2:
% Initiate IFE by calculating the factors and factor loading (beta=0)
if rfix==0
    [r2,~,~,~]=nbplog(Y,rmax,1,0);
elseif rfix==1
    r2=rmax;
end
[F12,L12,VNT2]=panelFactorNew(Y,r2);

%Perform IFE iteration
[betaiIFE(:,N+1:2*N), F2,L2, VNT2, e2, nnn(1,2),r2]=Mul_betaiIterNewET(X,XXinv, Y, F12,L12, r2,rmax,rfix, 0.00000001, betaIFEp, seIFEp);
sigma2(1,2)=trace(e2*e2')/(N*T-r2*(N+T)+r2^2-2);
SSR(1,2)=trace(e2*e2');
betaIFEH(:,2)=mean(betaiIFE0(:,N+1:2*N),2);
squaredIFEi=zeros(p,p,N);
for i=1:N
    squaredIFEi(:,:,i)=(betaiIFE0(:,N+i)-betaIFEH(:,2))*(betaiIFE0(:,N+i)-betaIFEH(:,2))';
end;
varcovIFEi=(1/(N*(N-1)))*sum(squaredIFEi,3);
seIFEi(:,2)=sqrt(diag(varcovIFEi));



betaiIFE;
seIFEi;
sigma2;



