%����Jacobi������Chebyshev���ٷ�
function [x,  iter, err, time] = mcheby( A, b, x, a1,b1,tol, max_it )
% ����: ϵ������A,�Ҷ�����b,��ʼ����x,�������tol, ���������� max_it
% �����������x,����err,��������iter,CPUʱ��time
if nargin<7, max_it=1000; end
if nargin<6, tol=1.e-6; end
if nargin<5, x=zeros(size(b)); end
iter = 0;   bnrm2 = norm(b);
if (bnrm2 == 0.0), bnrm2 = 1.0; end
r =b-A*x;  %�����ʼ�в�
err = norm(r) / bnrm2;
if ( err < tol ), return; end
n = length(b);
D=diag(diag(A)); B=D\(D-A);  f=D\b; %Jacobi��������
ga=2/(2-a1-b1); w1=(2-a1-b1)/(b1-a1); 
alpha=1.0/(4*w1*w1);  rho=2;  x0=x;
x1=ga*(B*x0+f)+(1-ga)*x0;   
tic
for iter = 1:max_it,    % ������ʼ
    x=rho*(ga*(B*x1+f)+(1-ga)*x1)+(1-rho)*x0;
    r = b-A*x;    %����в�
    err = norm(r) / bnrm2; 
    if ( err <= tol ), break; end
    x0=x1;  x1=x;
    rho=1.0/(1-alpha*rho);
end
time=toc;