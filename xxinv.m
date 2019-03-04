% Computes (X'X)^{-1} for a three dimensional matrix. This code was kindly provided to me by Bia via email correspondence.

function [xxinv] = XXinv(X) 
   [T,N,p]=size(X);
   xx=zeros(p,p);
   
   if p==1;
       xx(1,1)=trace(X'*X);
       
   end
   if p > 1
        for k=1:p;
            X1=X(:,:,k);
         
              for m=k:p;
                  X2=X(:,:,m);
                  xx(k,m)=trace(X1'*X2);
                     if k < m;
                          xx(m,k)=xx(k,m);
                     end
             end
        end
    end
  
   xxinv=inv(xx);