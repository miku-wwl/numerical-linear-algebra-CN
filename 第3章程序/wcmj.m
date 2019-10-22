%基于Richardson格式的整体校正法
function [x,iter,time,res]=wcmj(A,b,x,m,s,tol)
n=length(b); In=eye(n);   
D=diag(diag(A));
B=D\(D-A);  f=D\b; iter=0; 
res=norm(b-A*x)/norm(b);
X=zeros(n,m); R=zeros(n,m); 
X(:,1)=x; R(:,1)=res;
tic
while (iter<1000)
    iter=iter+1;
    for i=2:m
        X(:,i)=B*X(:,i-1)+f;
        R(:,i)=b-A*X(:,i);
    end
    for i=1:m
        if i~=s
            Qk(:,i)=R(:,i)-R(:,s);
        end
    end
    y=-Qk\R(:,s);
    alpha(s)=1-sum(y);
    for i=1:s-1
        alpha(i)=y(i);
    end
    for i=s:m-1
        alpha(i+1)=y(i);
    end
    x=zeros(n,1);
    for i=1:m
        x=x+alpha(i)*X(:,i);
    end
    res=norm(b-A*x)/norm(b);
    if res<tol, break; end
    X(:,1)=x;   R(:,1)=res;
end
time=toc;

