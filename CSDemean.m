% This function removes cross-sectionally demeans. I.e., removes time period fixed effects.

function [Xdot, Ydot] = CSDemean(X,Y,RemoveFixedEffect)
[T,N,p]=size(X);
y=Y;

if RemoveFixedEffect==1
    in=ones(1,N);
    ydotvec = zeros(T,N);
    a=mean(y,2);
    for t=1:T
        ydotvec(t,:) = y(t,:) - in*a(t,1);
    end
    
    X1=zeros(T,N,p);
    for j=1:p
        b=mean(X(:,:,j),2);
        for t=1:T
            X1(t,:,j) = X(t,:,j) - in*b(t,:);
        end
    end
    
    Ydot=ydotvec;
    Xdot=X1;
    
else
    
    Ydot=Y;
    Xdot=X;
end