% Performs bias correction, based on equations (23) and (24) in Bai (2009).

function [BETAbc]=biasife(Xdot,F1,L1,e1,beta);
[T,N,p]=size(Xdot);


MF=eye(T) - (F1*F1')/T;
a=L1*inv((L1'*L1)/N)*L1';


Xdotp=permute(Xdot,[1 3 2]);    

V1=zeros(T,p*N);
for i=1:N
    V1=V1 + ( (1/N)*kron(a(i,:),Xdotp(:,:,i)) );
end
V=zeros(T,p,N);
for i=1:N
    V(:,:,i)=V1(:,(p*(i-1)+1):(p*(i-1)+p));
end;


Zinner1=zeros(T,p*N);
for i=1:N
    Zinner1=Zinner1 + ( (1/N)*kron(a(i,:),MF*Xdotp(:,:,i)) );
end
Zinner=zeros(T,p,N);
for i=1:N
    Zinner(:,:,i)=Zinner1(:,(p*(i-1)+1):(p*(i-1)+p));
end;


Z=zeros(T,p,N);
for i=1:N
    Z(:,:,i)=MF*Xdotp(:,:,i) - Zinner(:,:,i);
end;


sigma2=e1.^2;
sigma2i=mean(sigma2,1);

innersumD0=zeros(p,p,N);
for i=1:N
        innersumD0(:,:,i)=Z(:,:,i)'*Z(:,:,i);
end;

D0=(1/(N*T))*sum(innersumD0,3);


innersumB=zeros(p,N);
for i=1:N
    innersumB(:,i)=(((Xdotp(:,:,i) - V(:,:,i))'*F1)/T)*inv((L1'*L1)/N)*L1(i,:)'*sigma2i(i);
end;

B=-inv(D0)*(1/N)*sum(innersumB,2);


sigma2t=mean(sigma2,2);

omega=diag(sigma2t);

innersumC=zeros(p,N);
for i=1:N
    innersumC(:,i)=(Xdotp(:,:,i)'*MF*omega*F1)*inv((L1'*L1)/N)*L1(i,:)';
end;

C=-inv(D0)*(1/(N*T))*sum(innersumC,2);

biascorrected=beta-(1/N)*B-(1/T)*C;

BETAbc=biascorrected;
