%Jacobi��������-mjacobi.m
function [lambda,Q]=Jacobi_eig(A,tol)
%��������Jacobi������ʵ�Գƾ���A��ȫ������ֵ����������.
%�������:AΪn�׶ԳƷ���, tolΪ�������.�������:lambdaΪ��
%��,�����ΪA������ֵ, QΪ����,��Ԫ��Ϊ����A��n����������
if nargin<2, tol=1e-6; end
[n]=size(A,1); Q=eye(n);
%����A�ķǶԽ�Ԫ����ֵ���Ԫ�����ڵ���p����q
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