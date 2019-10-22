function [Iter, Err, t]=ex612BJ(A,B) 
tic;
[n,n1]=size(A); I=eye(n);  
F=I;  X0=I; Iter=0; 
while(Iter<=1500)
    for i=1:n
        if(i==1)
            S2=zeros(n,1);
            for (j=i+1:n)
                S2=S2+B(j,i)*X0(:,j);
            end
            X(:,i)=(A+B(i,i)*I)\(F(:,i)-S2);
        elseif(i==n)
            S1=zeros(n,1);
            for (j=1:i-1)
                S1=S1+B(j,i)*X0(:,j);
            end
            X(:,i)=(A+B(i,i)*I)\(F(:,i)-S1);
        else
            S1=zeros(n,1); S2=zeros(n,1);
            for (j=1:i-1),   S1=S1+B(j,i)*X0(:,j);  end
            for (j=i+1:n),  S2=S2+B(j,i)*X0(:,j);  end
            X(:,i)=(A+B(i,i)*I)\(F(:,i)-S1-S2);
        end
    end
    Err=norm(X-X0,'fro');
    if(Err<=1.e-10)
       break;
    end
    X0=X;  Iter=Iter+1;
end
t=toc;