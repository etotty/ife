
% computing panel data estimator given the projection matrix M_F. This code was kindly provided to me by Bia via email correspondence.

function [beta] = Mul_panelbeta(X, Y, MF) 
   [T,N,p]=size(X);
   xx=zeros(p,p);
   xy=zeros(p,1);
   if p==1;
       xx(1,1)=trace(X'*MF*X);
       xy(1)=trace(X'*MF*Y);
   end
   if p > 1
        for k=1:p;
            MFX1=MF*X(:,:,k);
          %  xx(k,k)=trace(MFX1'*MFX1);
            xy(k)=trace(MFX1'*Y);
              for m=k:p;
                  MFX2=MF*X(:,:,m);
                  xx(k,m)=trace(MFX1'*MFX2);
                     if k < m;
                          xx(m,k)=xx(k,m);
                     end
             end
        end
    end
  
   beta=inv(xx)*xy;
    