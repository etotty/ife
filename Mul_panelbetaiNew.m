% Re-estimates heterogeneous OLS regression coefficients, given estimate of factor structure.


function [betai] = Mul_panelbetaiNew(X, xxinv, Y, F,L) 
   [T,N,p]=size(X);
   betai=zeros(p,N);

   for i=1:N
       xi=permute(X(:,i,:),[1 3 2]);
       yi=Y(:,i)-F*L(i,:)';
       betai(:,i)=inv(xi'*xi)*xi'*yi;
   end
    
