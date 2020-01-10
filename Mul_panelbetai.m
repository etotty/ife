

function [betai] = Mul_panelbetai(X, Y) 
   [T,N,p]=size(X);
   betai=zeros(p,N);

   for i=1:N
       i
       xi=permute(X(:,i,:),[1 3 2]);
       yi=Y(:,i);
       betai(:,i)=inv(xi'*xi)*xi'*yi;
   end
    