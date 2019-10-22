%Jacobi方法程序-mjacobi.m
function [lambda,Q]=Jacobi_eig(A,tol)
%本程序用Jacobi方法求实对称矩阵A的全部特征值和特征向量.
%输入参数:A为n阶对称方阵, tol为容许误差.输出参数:lambda为向
%量,其分量为A的特征值, Q为矩阵,其元素为矩阵A的n个特征向量
if nargin<2, tol=1e-6; end
[n]=size(A,1); Q=eye(n);
%计算A的非对角元绝对值最大元素所在的行p和列q
[w1,p]=max(abs(A-diag(diag(A))));
[w2,q]=max(w1); p=p(q);
while(1)
  d=(A(p,p)-A(q,q))/(2*A(p,q));
  if(d>=0)
     t=-d+sqrt(d^2+1);
  else
     t=-d-sqrt(d^2+1);
  end
  c=1/sqrt(t^2+1); s=c*t; G=[c s; -s c];
  A([p q],:)=G*A([p q],:);
  A(:,[p q])=A(:,[p q])*G';
  Q(:,[p q])=Q(:,[p q])*G';
  [w1,p]=max(abs(A-diag(diag(A))));
  [w2,q]=max(w1); p=p(q);
  if (abs(A(p,q))<tol*sqrt(sum(diag(A).^2)/n))
     break;
  end
end
lambda=sort(diag(A));