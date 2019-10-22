%Jacobi-Davidson方法程序-jacobidavidson.m
function [mu,u,Vm]=jacobidav(A,v,tol,max_it)
%本程序用Jacobi-Davidson求矩阵A的模最大的特征值及相应的特征向量
%输入参数:A为n阶实方阵,v为初始向量,tol为容许误差,max_it为子空间的最大维数.
%输出参数:mu返回A的模最大特征值,u为相应的特征向量,Vm为子空间基矩阵,满足Vm'*Vm=Im
if nargin<4, max_it=ceil(size(A,1)/2); end
if nargin<3, tol=1e-6; end
n=size(A,1);  m=0; Vm=[ ];  I=eye(n);
v=v/norm(v);
Vm=[Vm,v]; Am=Vm'*A*Vm;
while(m<max_it)
    m=m+1;
    [Y,D]=eig(Am); [s,t]=max(abs(diag(D)));
    mu=D(t,t); u=Vm*Y(:,t); r=A*u-mu*u;
    if (norm(r)<=tol), break; end
    At=(I-u*u')*(A-mu*I)*(I-u*u');
    [z]=gmres(At,-r); %系统自带的GMRES函数
    v=z+u; v=v/norm(v); 
    X=[Vm,v]; [Vm]=House_orth(X);  
    %Am=[Am, Vm'*A*v; v'*A*Vm, v'*A*v];
    %X=Vm; [Vm]=G_Schmidt2(X);
    Am=Vm'*A*Vm;
    norm(Vm'*Vm-eye(m+1))
end
m

   