%����-mddqrdec.m %˫��λ����ʽQR����
function [iter,D]=ddqrdec(A,tol)
%��;: ��˫��λ����ʽQR������ʵ�����ȫ������ֵ.
%����: n��ʵ����A, ���ƾ���tol(Ĭ����1.e-5)
%���: ��������Iter, A��ȫ������ֵD
if nargin<2,tol=1e-5;end
n=size(A,1);
D=zeros(n,1); i=n; m=n; iter=0;  %��ʼ��
A=mhessen(A);  %������AΪHessenberg����
while (m>0)   
    %��˫��λ����ʽQR�������е���
   if m<=2
      la=eig(A(1:m,1:m)); 
      D(1:m)=la';  break;
   end
   iter=iter+1;
   A=ddqrtran(A,m); 
   %����Hessenberg ������QR�ֽ�,�����������Ʊ任
   %����ĳ�������ж��Ƿ���ֹ
   for k=m-1:-1:1
      if abs(A(k+1,k))<tol
          if m-k<=2
             la=eig(A(k+1:m,k+1:m));
             j=i-m+k+1; D(j:i)=la';
             i=j-1; m=k; break;
           end
       end
   end
end
