%求矩阵全部特征值的基本QR方法程序-qr_eig.m
function [iter,D]=qr_eig(A,tol,N)
%本程序用基本QR算法求n阶实方阵A的全部特征值,其中tol为
%控制精度,N为最大迭代次数,输出迭代次数Iter和A的全部特征值D
%调用函数: mhessen.m,hessen_qrtran.m,eig-仅用于1,2矩阵
if nargin<3,N=1000;end
if nargin<2,tol=1e-5;end
n=size(A,1); D=zeros(n,1); 
i=n; m=n; iter=0;  %初始化
A=mhessen(A);  %化矩阵A为Hessenberg矩阵
while (iter<=N)   %用基本QR算法进行迭代
    if m<=2
        la=eig(A(1:m,1:m)); D(1:m)=la';  
        break;
    end
    iter=iter+1;
   %对上Hessenberg矩阵做QR分解并做正交相似变换
   A=hessen_qrtran(A,m); 
   %下面的程序段是判断是否终止
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
