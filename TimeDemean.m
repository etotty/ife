% This program time demeans the data. I.e., removes unit fixed effects.

function [Xdot, Ydot] = TimeDemean(X,Y,RemoveFixedEffect)
[T,N,p]=size(X);
y=Y;

if RemoveFixedEffect==1
    it=ones(T,1);
    ydotvec = zeros(T,N);
    a=mean(y);
    for i=1:N
        ydotvec(:,i) = y(:,i) - it*a(i);
    end
    
    X1=zeros(T,N,p);
    for j=1:p
        b=mean(X(:,:,j));
        for i=1:N
            X1(:,i,j) = X(:,i,j) - it*b(i);
        end
    end
    
    Ydot=ydotvec;
    Xdot=X1;
    
else
    
    Ydot=Y;
    Xdot=X;
end