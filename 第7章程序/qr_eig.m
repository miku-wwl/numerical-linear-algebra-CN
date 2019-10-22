%�����ȫ������ֵ�Ļ���QR��������-qr_eig.m
function [iter,D]=qr_eig(A,tol,N)
%�������û���QR�㷨��n��ʵ����A��ȫ������ֵ,����tolΪ
%���ƾ���,NΪ����������,�����������Iter��A��ȫ������ֵD
%���ú���: mhessen.m,hessen_qrtran.m,eig-������1,2����
if nargin<3,N=1000;end
if nargin<2,tol=1e-5;end
n=size(A,1); D=zeros(n,1); 
i=n; m=n; iter=0;  %��ʼ��
A=mhessen(A);  %������AΪHessenberg����
while (iter<=N)   %�û���QR�㷨���е���
    if m<=2
        la=eig(A(1:m,1:m)); D(1:m)=la';  
        break;
    end
    iter=iter+1;
   %����Hessenberg������QR�ֽⲢ���������Ʊ任
   A=hessen_qrtran(A,m); 
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
