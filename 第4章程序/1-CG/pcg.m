%Ԥ�������ݶȷ�����-pcg.m
function [x,iter,time,res,resvec]=pcg(A,b,x,M,max_it,tol)
%����:ϵ������A,�Ҷ�����b,��ʼ����x,Ԥ��
%����M,�������tol����������max_it
%���:������x,��������iter,CPUʱ��time,
%��ֹʱ��Բв�ģres��Բв�ģ����resvec
tic; r=b-A*x; z=M\r; p=z; 
rho=z'*r; mr=norm(r); iter=0;
while (iter<max_it)
    iter=iter+1;
    u=A*p; alpha=rho/(p'*u);
    x=x+alpha*p; r=r-alpha*u;
    z=M\r; rho1=z'*r;  
    beta=rho1/rho; p=z+beta*p;
    res=norm(r)/mr;
    resvec(iter)=res;
    if (res<tol), break; end
    rho=rho1;
end
time=toc;
