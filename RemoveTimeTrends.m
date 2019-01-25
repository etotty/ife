%This function removes unit-level or S-group-level time trends.

function [Xdetrend, Ydetrend] = RemoveTimeTrends(X1, Y1, order, S)
[T,N,p]=size(X1);
m=length(S(:,1));
trend=[1:T]';

if order>0
    
    if m==1; %removes time trend for each of the "N" cross-section units
        Xi=zeros(T,order+1);
        Xi(:,1)=ones(T,1);
        if order>0
            for o=1:order
                Xi(:,o+1)=[trend.^o];
            end
            betatrend=zeros(order+1,N);
            Ydetrend=zeros(T,N);
            for i=1:N
                Yi=Y1(:,i);
                betatrend(:,i)=inv(Xi'*Xi)*Xi'*Yi;
                Ydetrend(:,i)=Yi-Xi*betatrend(:,i);
            end
        else
            Ydetrend=Y1;
        end
        
        Xdetrend=zeros(T,N,p);
        for j=1:p
            if order>0
                for o=1:order
                    Xi(:,o+1)=[trend.^o];
                end
                betatrend=zeros(order+1,N);
                for i=1:N
                    x1i=X1(:,i,j);
                    betatrend(:,i)=inv(Xi'*Xi)*Xi'*x1i;
                    Xdetrend(:,i,j)=x1i-Xi*betatrend(:,i);
                end
            else
                Xdetrend(:,:,j)=X1(:,:,j);
            end
        end
        
        
    else %places each of N cross-section units into one of "m" groups (m<N) and removes time trend for each m
        
        m=length(S(1,:));
        Ydetrend=zeros(T,N);
        Xdetrend=zeros(T,N,p);
        for s=1:m
            nn=sum(S(:,s));
            Xs=zeros(T,nn,p);
            Ss=zeros(N,nn);
            count=0;
            for i=1:N
                count=count+S(i,s);
                if S(i,s)>0
                    Ss(i,count)=1;
                end
            end
            Ys=Y1*Ss;
            for j=1:p
                Xs(:,:,j)=X1(:,:,j)*Ss;
            end
            
            XX=zeros(order+1,order+1,nn);
            XY=zeros(order+1,1,nn);
            for i=1:nn
                Xi=zeros(T,order+1);
                Xi(:,1)=ones(T,1);
                for o=1:order
                    Xi(:,o+1)=[trend.^o];
                end
                XX(:,:,i)=(1/nn)*(Xi'*Xi);
                Yi=Ys(:,i);
                XY(:,:,i)=(1/nn)*(Xi'*Yi);
            end
            XXsum=sum(XX,3);
            XYsum=sum(XY,3);
            betatrend=inv(XXsum)*XYsum;
            
            Ysdetrend=zeros(T,nn);
            for i=1:nn
                Ysdetrend(:,i)=Ys(:,i)-Xi*betatrend
            end
            
            Ydetrend=Ydetrend+Ysdetrend*Ss';
            
            Xsdetrend=zeros(T,nn,p);
            for j=1:p
                XX=zeros(order+1,order+1,nn);
                Xx=zeros(order+1,1,nn);
                for i=1:nn
                    Xi=zeros(T,order+1);
                    Xi(:,1)=ones(T,1);
                    for o=1:order
                        Xi(:,o+1)=[trend.^o];
                    end
                    XX(:,:,i)=(1/nn)*(Xi'*Xi);
                    xi=Xs(:,i,j);
                    Xx(:,:,i)=(1/nn)*(Xi'*xi);
                end
                XXsum=sum(XX,3);
                Xxsum=sum(Xx,3);
                betatrend=inv(XXsum)*Xxsum;
                
                for i=1:nn
                    Xsdetrend(:,i,j)=Xs(:,i,j)-Xi*betatrend
                end
                
                Xdetrend(:,:,j)=Xdetrend(:,:,j)+Xsdetrend(:,:,j)*Ss';
            end
            
            
            
        end
    end
else
    Ydetrend=Y1;
    Xdetrend=X1;
end























