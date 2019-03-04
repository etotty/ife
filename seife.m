% Creates se for IFE using variables as defined in Bai (2009), one variable at a time (p=1)

function [se]=seIFE(Xdot,F1,L1,e1);
[T,N,p]=size(Xdot);

MF=eye(T) - (F1*F1')/T;

a=L1*inv((L1'*L1)/N)*L1';

Xdotp=permute(Xdot,[1 3 2]);

inner=zeros(T,p*N);
for k=1:N
    inner=inner + ( (1/N)*kron(a(k,:),MF*Xdotp(:,:,k)) );
end

Z=zeros(T,p,N);
for i=1:N
    Z(:,:,i)=MF*Xdotp(:,:,i) - inner(:,(p*(i-1)+1):(p*(i-1)+p));
end
    
innersumD0=zeros(p,p,N);
for i=1:N
        innersumD0(:,:,i)=Z(:,:,i)'*Z(:,:,i);
end

D0itsum=sum(innersumD0,3);
D0=(1/(N*T))*D0itsum;

innersumD3=zeros(p,p,T,N);
sigma2=e1.^2;
for i=1:N
    for t=1:T
        innersumD3(:,:,t,i)=Z(t,:,i)'*Z(t,:,i)*sigma2(t,i);
    end
end

D3isum=sum(innersumD3,4);
D3itsum=sum(D3isum,3);
D3=(1/(N*T))*D3itsum;

var=(1/(N*T))*inv(D0)*D3*inv(D0);
se=sqrt(diag(var));

