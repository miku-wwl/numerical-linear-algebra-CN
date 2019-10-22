%Jacobi-Davidson方法程序--jacobidavidson.m
function [mu,u,Vm]=mjacobidav(A,v1,tol,max_it)
%本程序用Jacobi-Davidson求矩阵A的某特征值特征值,max_it为子空间的最大维数.
%输入参数:A为n阶实对称方阵, tol为容许误差.
%输出参数:mu返回A的某个特征值
if nargin<4, max_it=size(A,1); end
if nargin<5, tol=1e-5; end
n=size(A,1);  m=0; Vm=[ ]; I=eye(n);
Vm=[Vm,v1];   Am=Vm'*A*Vm;
while(m<max_it)
    m=m+1;   
    if m<=2
        [Y,d]=eig(Am);
        [ss,k]=max(abs(diag(d)));
        mu=d(k,k); u=Vm*Y(:,k);
    else
        x=ones(size(Am,1),1);
        [mu,y,k1]=mypower(Am,x);
        [mu,y,k2]=mvpower(Am,x,mu);
        u=Vm*y;
    end
    
    r=A*u-mu*u;
    if (norm(r)<=tol), break; end
    At=(I-u*u')*(A-mu*I)*(I-u*u');
    rt=(A-mu*I)*u;
    z=-At\rt;   v1=z;
    for i=1:m
        v1=v1-(z'*Vm(:,i))*Vm(:,i);
    end
    v1=v1/norm(v1);
    Am=[Am, Vm'*A*v1; v1'*A*Vm, v1'*A*v1];
    %pause
    Vm=[Vm,v1];
end
   