%Jacobi-Davidson方法程序--jacobidavidson.m
function [mu,u,Vm]=jacobidavi(A,v,tol,max_it)
%本程序用Jacobi-Davidson求矩阵A的某特征值特征值
%输入参数:A为n阶实方阵,v1为初始向量,
%tol为容许误差,max_it为子空间的最大维数.
%输出参数:mu返回A的某个特征值,u为相应的
%特征向量,Vm为子空间基矩阵,满足Vm'*Vm=Im
if nargin<4, max_it=ceil(size(A,1)/2); end
if nargin<3, tol=1e-5; end
n=size(A,1);  m=0; Vm=[ ];  I=eye(n);
x0=rand(n,1); v=v/norm(v);
Vm=[Vm,v];   Am=Vm'*A*Vm;
while(m<max_it)
    m=m+1;
    [Y,Lam]=eig(Am);
    [ss,t]=max(abs(diag(Lam)));
    mu=Lam(t,t);
    u=Vm*Y(:,t); r=A*u-mu*u;
    if (norm(r)<=tol), break; end
    %At=A-mu*I; 
    %p1=gmres(At, r); p2=gmres(At, u); %p1=At\r; p2=At\u;
    %alpha=(u'*p1)/(u'*p2); f=-r+alpha*u;
    %z=(A-mu*I)\(-r+alpha*u);
    %[z] = gmres(At, f);  
    At=(I-u*u')*(A-mu*I)*(I-u*u');
    %M=diag(diag(At));
    %[x,flag,relres,iter,resvec] = gmres(A,b,restart,tol,maxit,M1,M2,x0)
    [z] = gmres(At,-r);   %系统自带的GMRES函数
    %z=-At\r; 
    s=zeros(n,1); v=z;
    for i=1:m
        v=v-(z'*Vm(:,i))*Vm(:,i);
    end
    %v=z-s; 
    v=v/norm(v);
    %Am=[Am, Vm'*A*v; v'*A*Vm, v'*A*v];
    %pause
    Vm=[Vm,v]; 
%     X=[Vm,v];
%     [Vm]=G_Schmidt2(X);
%     [Vm]=G_Schmidt2(Vm);
    Am=Vm'*A*Vm;
    norm(Vm'*Vm-eye(m+1))
end

   