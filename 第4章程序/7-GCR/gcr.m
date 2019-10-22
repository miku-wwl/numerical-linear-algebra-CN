%广义共轭残量法(GCR方法)程序-gcr.m
function [x,k,time,res,resvec]=gcr(A,b,x,max_it,tol)
tic; n=length(b); r=b-A*x; p=r; mr=norm(r); 
resvec(1)=1; AP=zeros(n,max_it+1); 
P=zeros(n,max_it+1); 
AP(:,1)=A*p; P(:,1)=p;  k=0;
while (k<max_it)
    k=k+1;
    alpha=(r'*AP(:,k))/(AP(:,k)'*AP(:,k));
    x=x+alpha*P(:,k);
    r=r-alpha*AP(:,k);
    Ar=A*r; s1=zeros(n,1); s2=zeros(n,1);
    for (i=1:k)
        b(i)=-(Ar'*AP(:,i))/(AP(:,i)'*AP(:,i));
        s1=s1+b(i)*P(:,i);
        s2=s2+b(i)*AP(:,i);
    end
    P(:,k+1)=r+s1;
    AP(:,k+1)=Ar+s2;
    res=norm(r)/mr; resvec(k+1)=res;
    if (res<tol),   break;   end
end
time=toc;
