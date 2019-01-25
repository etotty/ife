% Re-estimates regression coefficients, given estimate of factor structure. This code was kindly provided to me by Bia via email correspondence.


function [beta] = Mul_panelbetaNew(X, xxinv, Y, F,L) 
   [T,N,p]=size(X);
  
   xy=zeros(p,1);
  
   for k=1:p;
            
          
   xy(k)=trace(  X(:,:,k)'*(Y-F*L')  );
              
    end
  
   beta=xxinv*xy;
    