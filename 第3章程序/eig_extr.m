%����Jacobi����������ֵ���Ʒ�
function [x,iter,res,time]=eig_extr(A,b,x,m,lamd,tol)
n=length(b); In=eye(n);   
D=diag(diag(A)); B=D\(D-A);  f=D\b;  %Jacobi��������
iter=0; 
res=norm(b-A*x)/norm(b);
X=zeros(n,m); R=zeros(n,m); 
X(:,1)=x; R(:,1)=res;
tic
while (iter<1000)
    iter=iter+1;
    for k=2:m
        X(:,k)=B*X(:,k-1)+f;
        R(:,k)=b-A*X(:,k);
    end
    x=X(:,m-1)+(X(:,m)-X(:,m-1))/(1-lamd);
    res=norm(b-A*x)/norm(b);
    if res<tol, break; end
    X(:,1)=x;   R(:,1)=res;
end
time=toc;

