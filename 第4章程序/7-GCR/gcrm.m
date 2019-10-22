%GCR(m)·½·¨³ÌÐò-gcrm.m
function [x,out,int,time,res,resvec]=gcrm(A,b,x,restrt,max_it,tol)
tic; n=length(b); m=restrt;
r=b-A*x; p=r; mr=norm(r); k=0;
AP=zeros(n,m+1); P=zeros(n,m+1);
AP(:,1)=A*p; P(:,1)=p; resvec(1)=1;
while (k<max_it)
    k=k+1;
    for i=1:m
        alpha=(r'*AP(:,i))/(AP(:,i)'*AP(:,i));
        x=x+alpha*P(:,i);
        r=r-alpha*AP(:,i);
        res=norm(r)/mr;        
        resvec((k-1)*m+1+i)=res;
        if (res<tol), break;  end
        Ar=A*r; s1=zeros(n,1); s2=zeros(n,1);
        for (s=1:i)
            b(s)=-(Ar'*AP(:,s))/(AP(:,s)'*AP(:,s));
            s1=s1+b(s)*P(:,s);
            s2=s2+b(s)*AP(:,s);
        end
        P(:,i+1)=r+s1;
        AP(:,i+1)=Ar+s2;
    end
    if (res<tol), 
         out=k; int=i;  break;   
    end
    AP(:,1)=AP(:,m+1); P(:,1)=P(:,m+1);
end
time=toc;
