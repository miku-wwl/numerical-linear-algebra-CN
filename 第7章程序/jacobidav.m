%Jacobi-Davidson��������-jacobidavidson.m
function [mu,u,Vm]=jacobidav(A,v,tol,max_it)
%��������Jacobi-Davidson�����A��ģ��������ֵ����Ӧ����������
%�������:AΪn��ʵ����,vΪ��ʼ����,tolΪ�������,max_itΪ�ӿռ�����ά��.
%�������:mu����A��ģ�������ֵ,uΪ��Ӧ����������,VmΪ�ӿռ������,����Vm'*Vm=Im
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
    [z]=gmres(At,-r); %ϵͳ�Դ���GMRES����
    v=z+u; v=v/norm(v); 
    X=[Vm,v]; [Vm]=House_orth(X);  
    %Am=[Am, Vm'*A*v; v'*A*Vm, v'*A*v];
    %X=Vm; [Vm]=G_Schmidt2(X);
    Am=Vm'*A*Vm;
    norm(Vm'*Vm-eye(m+1))
end
m

   