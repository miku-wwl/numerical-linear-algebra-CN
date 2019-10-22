function [Iter, Err, t]=ex612GS(A,B) 
tic; 
DA=diag(diag(A));  [n,n1]=size(DA);
I=eye(n);  F=I;  X0=I; Iter=0; 
while(Iter<=1500)
    for i=1:n
        if (i==1)
            for(l=1:n) 
                if(l==1)
                    X(i,l)=(F(i,l) - A(l,l+1:n)*X0(i,l+1:n)'- B(i+1:n,i)'*X0(i+1:n,l))/(B(i,i)+A(l,l));                    
                elseif(l==n)
                    X(i,l)=(F(i,l) - A(l,1:l-1)*X(i,1:l-1)' - B(i+1:n,i)'*X0(i+1:n,l))/(B(i,i)+A(l,l));        
                else
                    X(i,l)=(F(i,l) - A(l,1:l-1)*X(i,1:l-1)' - A(l,l+1:n)*X0(i,l+1:n)' - B(i+1:n,l)'*X0(i+1:n,l))/(B(i,i)+A(l,l));
                end
            end
        elseif(i==n)
            for(l=1:n) 
                if(l==1)
                    X(i,l)=(F(i,l) - B(1:i-1,i)'*X(1:i-1,l) - A(l,l+1:n)*X0(i,l+1:n)'-B(i+1:n,i)'*X0(i+1:n,l))/(B(i,i)+A(l,l)); 
                elseif(l==n)
                    X(i,l)=(F(i,l) - B(1:i-1,i)'*X(1:i-1,l) - A(l,1:i-1)*X(i,1:i-1)')/(B(i,i)+A(l,l)); 
                else
                    X(i,l)=(F(i,l) - B(1:i-1,i)'*X(1:i-1,l) - A(l,1:i-1)*X(i,1:i-1)' - A(l,l+1:n)*X0(i,l+1:n)')/(B(i,i)+A(l,l));
                end
            end
        else
            for(l=1:n)
                if(l==1)
                    X(i,l)=(F(i,l) - B(1:i-1,i)'*X(1:i-1,l) - A(l,l+1:n)*X0(i,l+1:n)' - B(i+1:n,i)'*X0(i+1:n,l))/(B(i,i)+A(l,l));   
                elseif(l==n)
                    X(i,l)=(F(i,l) - B(1:i-1,i)'*X(1:i-1,l) - A(l,1:i-1)*X(i,1:i-1)' - B(i+1:n,i)'*X0(i+1:n,l))/(B(i,i)+A(l,l));   
                else
                    X(i,l)=(F(i,l) - B(1:i-1,i)'*X(1:i-1,l) - A(l,1:i-1)*X(i,1:i-1)' - A(l,l+1:n)*X0(i,l+1:n)' - B(i+1:n,i)'*X0(i+1:n,l) )/(B(i,i)+A(l,l));
                end
            end
        end
    end
    Err=norm(X-X0,'fro');
    if(Err<=1.e-10)
       break;
    end
    X0=X;  Iter=Iter+1;
end
t=toc;