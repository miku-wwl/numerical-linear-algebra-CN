%重开始PGMRES方法-PGMRES(m)-gmresm.m
function [x,out,int,time,res,resvec,flag]=pgmresm(A,b,x,M,restrt,max_it,tol )
tic; flag=0; int=0; r=M\(b-A*x);  %计算残差
beta=norm(r); res=norm(r)/beta; resvec(1)=res;
n =length(b);  m=restrt;
V(1:n,1:m+1)=zeros(n,m+1);
H(1:m+1,1:m)=zeros(m+1,m);
c(1:m)=zeros(m,1); s(1:m)=zeros(m,1);
e1=zeros(n,1); e1(1)=1.0;
for k=1:max_it,  
    V(:,1)=r/norm(r); xi=norm(r)*e1;
    for j=1:m,  %用Arnoldi方法构造正交基 
        w=M\(A*V(:,j));     
        for i =1:j
            H(i,j)=w'*V(:,i);  w=w-H(i,j)*V(:,i);
        end
        H(j+1,j)=norm(w);
        if abs(H(j+1,j))/beta<tol,
            return;
        else
            V(:,j+1)=w/H(j+1,j);
        end
        for i=1:j-1, %第i次Givens变换
            temp=c(i)*H(i,j)+s(i)*H(i+1,j);
            H(i+1,j)=-s(i)*H(i,j)+c(i)*H(i+1,j);
            H(i,j)=temp;
        end
        [c(j),s(j),H(j,j)]=givens(H(j,j),H(j+1,j) ); %第j次Givens变换
        xi(j+1)=-s(j)*xi(j); xi(j)=c(j)*xi(j); H(j+1,j)=0.0;
        res=abs(xi(j+1))/beta; resvec((k-1)*m+j+1)=res;
        if (res<=tol),                         
            y=H(1:j,1:j)\xi(1:j); x=x+V(:,1:j)*y; 
            break; %跳出内循环
        end
    end
    if (res<tol )
        out=k; int=j; break;  %跳出外循环
    end
    y=H(1:m,1:m)\xi(1:m);
    x=x+V(:,1:m)*y;                          
    r=M\(b-A*x);       
end
if (res>tol), flag=1; end;  %不收敛
time=toc;