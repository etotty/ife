%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produces Bai's (2009) interactive fixed effects (IFE) estimator (equations 11-12), including 
% bias-correction for correlations and heteroskedasticities (section 7), and 
% standard errors (section 7-8).
% Bai, Jushan. 2009. "Panel Data Models With Interactive Fixed Effects." Econometrica, 77 (4), 1229-1279.
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

% OUTPUTS:
% betaIFE = (px2) vector of IFE regression coefficients for each starting method
% seIFE = (px2) vector of standard errors for IFE regression coefficients for each starting methods
% sigma2 = (1,2) vector of estimated residual variance for each starting method
% SSR = (1,2) vector of final sum of squared residuals for each starting method
% nnn = (1,2) vector of number of iterations before convergence for each starting method
% r1 = final estimate of number of factors for starting method 1
% r2 = final estimate of number of factors for starting method 2
% betaIFE0 = (px1) vector of IFE regression coefficients before performing bias correction step
% betaOLS = (px1) vector of OLS regression coefficients
% F1 = estimated factors for starting method 1
% L1 = estimated loadings for starting method 1
% F2 = estimated factors for starting method 2
% L2 = estiamted loadings for starting method 2
% VNT1 = diagnonal matrix with the r-largest eigenvalues of the second moment matrix of the regression residuals
% VNT2 = diagnonal matrix with the r-largest eigenvalues of the second moment matrix of the regression residuals



function [betaIFE, seIFE, sigma2, SSR, nnn, r1, r2, betaIFE0, betaOLS, F1, L1, F2, L2, VNT1, VNT2] = IFEfunction(X1,Y1,order,S,demean,rmax,rfix); 

[T,N,p]=size(X1);


%%%%% 1. Remove state-specific trends and demean
[Xdet,Ydet] = RemoveTimeTrends(X1,Y1,order,S);
% remove fixed effects
if demean==0
    Xdot=Xdet;
    Ydot=Ydet;
elseif demean==1
    [Xdot,Ydot]=CSDemean(Xdet,Ydet,1);
elseif demean==2
    [Xdot,Ydot]=TimeDemean(Xdet,Ydet,1);
elseif demean==3
    [Xdot,Ydot]=DemeanData(Xdet,Ydet);
end
Y=Ydot;
X=Xdot;


%%%%% 2. Estimation
% calculate XXinv outside of loop to save time
[XXinv]= xxinv(X); 

betaIFE0=zeros(p,2);
betaIFE=zeros(p,2);    % to contain the interative effect estimators, with different staring methods (2 methods)
nnn=zeros(1,2);   % contain the number of iterations for the 2 methods to achieve convergence
sigma2=zeros(1,2); % the estimated residual variance for each starting method, also the optimal value of the objective function
seIFE=zeros(p,2);     % contains the IFE standard errors for the 2 different starting methods
SSR=zeros(1,2);

%%%%% 2a. Starting method 1:
% Initiate IFE by calculating standard panel beta estimator (lambda=0) as a starting value
[betaOLS]=Mul_panelbeta(X,Y,eye(T)); 

U=Y;
for k=1:p;
    U=U-X(:,:,k)* betaOLS(k,1);
end
U2=U.^2;
sigsq=trace(U*U')/(N*T-p);
varcov=sigsq*XXinv*inv(XXinv)*XXinv;
seOLS=sqrt(diag(varcov));

% Get the first estimate of r and F,lambda (because iteration below starts by estimating beta
if rfix==0
    [r1,~,~,~]=nbplog(U,rmax,1,0);
elseif rfix==1
    r1=rmax;
end
[F1,L1,VNT]=panelFactorNew(U,r1);

% Perform IFE iteration
[betaIFE0(:,1), F1,L1, VNT1, e1, nnn(1,1),r1]=Mul_betaIterNew(X,XXinv, Y, F1,L1, r1,rmax,rfix, 0.00000001);
sigma2(1,1)=trace(e1*e1')/(N*T-r1*(N+T)+r1^2-2);
SSR(1,1)=trace(e1*e1');
seIFE(:,1)=seife(X,F1,L1,e1);

% correct beta for serial corr, cross-sec corr, and heteroskedasticity
betaIFE(:,1)=biasife(X,F1,L1,e1,betaIFE0(:,1));


%%%%% 2b. Starting method 2:
% Initiate IFE by calculating the factors and factor loading (beta=0)
if rfix==0
    [r2,~,~,~]=nbplog(Y,rmax,1,0);
elseif rfix==1
    r2=rmax;
end
[F12,L12,VNT2]=panelFactorNew(Y,r2);

%Perform IFE iteration
[betaIFE0(:,2), F2,L2, VNT2, e2, nnn(1,2),r2]=Mul_betaIterNew(X,XXinv, Y, F12,L12, r2,rmax,rfix, 0.00000001);
sigma2(1,2)=trace(e2*e2')/(N*T-r2*(N+T)+r2^2-2);
SSR(1,2)=trace(e2*e2');
seIFE(:,2)=seife(X,F2,L2,e2);

% correct beta for serial corr, cross-sec corr, and heteroskedasticity
betaIFE(:,2)=biasife(X,F2,L2,e2,betaIFE0(:,2));




betaIFE;
seIFE;
sigma2;


