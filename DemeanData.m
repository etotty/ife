% This program demeans the variables in both directions. I.e., removes time and unit fixed effects.

function [Xdot, Ydot] = DemeanData(X,Y)
[T,N,p]=size(X);
y=Y;

ydotvec = zeros(T,N);
a=mean(y);
b=mean(y,2);
c=mean(y(:));
it = ones(T,1);
for i=1:N
    ydotvec(:,i) = y(:,i) - it*a(i) - b + it*c;
end

X1=zeros(T,N,p);
for j=1:p
    a=mean(X(:,:,j));
    b=mean(X(:,:,j),2);
    xx=X(:,:,j);
    c=mean(xx(:));
    it = ones(T,1);
    for i=1:N
        X1(:,i,j) = X(:,i,j) - it*a(i) - b + it*c;
    end
end



Ydot=ydotvec;
Xdot=X1;