function [iter,D]=ddiqr_eig(A,tol)
%��������˫��λ����ʽQR������ʵ�����ȫ������ֵ.
%����n��ʵ����A, ���ƾ���tol(Ĭ����1e-5)
%�����������iter, A��ȫ������ֵD
if nargin<2,tol=1e-5;end
n=size(A,1);
D=zeros(n,1); i=n; m=n; iter=0;  %��ʼ��
[A]=hessenb(A); %������AΪHessenberg����
while (m>0) 
    %��˫��λ����ʽQR�������е���
    if m<=2
        la=eig(A(1:m,1:m)); 
        D(1:m)=la';  break;
    end
    iter=iter+1;
    A=ddiqr_tran(A,m);  %����Hessenberg ������QR�ֽ�,�����������Ʊ任
    for k=m-1:-1:1  %����ĳ�������ж��Ƿ���ֹ
        if abs(A(k+1,k))<tol
            if m-k<=2
                la=eig(A(k+1:m,k+1:m));
                j=i-m+k+1; D(j:i)=la';
                i=j-1; m=k; break;
            end
        end
    end
end
