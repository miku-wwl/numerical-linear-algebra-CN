%基于Jacobi迭代的Chebyshev加速法
function [x,  iter, err, time] = mcheby( A, b, x, a1,b1,tol, max_it )
% 输入: 系数矩阵A,右端向量b,初始向量x,容许误差tol, 最大迭代次数 max_it
% 输出：解向量x,误差范数err,迭代次数iter,CPU时间time
if nargin<7, max_it=1000; end
if nargin<6, tol=1.e-6; end
if nargin<5, x=zeros(size(b)); end
iter = 0;   bnrm2 = norm(b);
if (bnrm2 == 0.0), bnrm2 = 1.0; end
r =b-A*x;  %计算初始残差
err = norm(r) / bnrm2;
if ( err < tol ), return; end
n = length(b);
D=diag(diag(A)); B=D\(D-A);  f=D\b; %Jacobi迭代矩阵
ga=2/(2-a1-b1); w1=(2-a1-b1)/(b1-a1); 
alpha=1.0/(4*w1*w1);  rho=2;  x0=x;
x1=ga*(B*x0+f)+(1-ga)*x0;   
tic
for iter = 1:max_it,    % 迭代开始
    x=rho*(ga*(B*x1+f)+(1-ga)*x1)+(1-rho)*x0;
    r = b-A*x;    %计算残差
    err = norm(r) / bnrm2; 
    if ( err <= tol ), break; end
    x0=x1;  x1=x;
    rho=1.0/(1-alpha*rho);
end
time=toc;