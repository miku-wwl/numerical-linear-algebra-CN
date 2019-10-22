%PGMRES方法程序-pgmres.m
function [x,k,time,res,resvec,flag]=pgmres(A,b,x,M,max_it,tol )
tic; flag=0; r=M\(b-A*x); beta=norm(r);  %计算残差
n=length(b); e1=zeros(n,1); e1(1)=1.0;
res=norm(r)/beta; resvec(1)=res;
V(:,1)=r/beta; xi=beta*e1; k=0;
while (k<=max_it)   
    k=k+1;
    w=M\(A*V(:,k));   
    for i=1:k %修正Arnoldi过程
        H(i,k)=w'*V(:,i);
        w=w-H(i,k)*V(:,i);
    end
    H(k+1,k)=norm(w);
    if abs(H(k+1,k))/beta<tol,
        return;
    else
        V(:,k+1)=w/H(k+1,k);
    end
    for i=1:k-1,                            
        temp=c(i)*H(i,k) + s(i)*H(i+1,k);
        H(i+1,k)=-s(i)*H(i,k) +c(i)*H(i+1,k);
        H(i,k)=temp;
    end
    [c(k),s(k),H(k,k)]=givens(H(k,k), H(k+1,k) ); %第k次Givens变换
    xi(k+1)=-s(k)*xi(k); xi(k)=c(k)*xi(k); H(k+1,k) = 0.0;
    res=abs(xi(k+1))/beta; resvec(k+1)=res;
    if (res<=tol ),   
        y=H(1:k,1:k)\xi(1:k); x=x+V(:,1:k)*y; 
        break; %跳出循环
    end
end
if (res>tol), flag=1; end; %不收敛
time=toc;