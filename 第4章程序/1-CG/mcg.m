%�����ݶȷ�����-mcg.m
function [x,iter,time,res,resvec]=mcg(A,b,x,max_it,tol)
%����:ϵ����A,�Ҷ���b,��ʼֵx,�������tol,��������max_it
%���:������x,��������iter,CPUʱ��time,
%��ֹʱ��Բв�ģres��Բв�ģ����resvec
tic; r=b-A*x; p=r; 
rho=r'*r; mr=sqrt(rho); iter=0;
while (iter<max_it)
    iter=iter+1;
    z=A*p; alpha=rho/(z'*p);
    x=x+alpha*p; r=r-alpha*z;
    rho1=r'*r; beta=rho1/rho;
    p=r+beta*p;
    res=sqrt(rho1)/mr;
    resvec(iter)=res;
    if (res<tol), break; end
    rho=rho1;
end
time=toc;
