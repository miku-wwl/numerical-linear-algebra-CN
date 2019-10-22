function A=ddiqr_tran(A,n)
%双重位移隐式QR变换. 给定不可约n阶上Hessenberg矩阵A, 其最后2x2阶
%主子阵有特征值a和b. 本算法计算Q^TAQ并覆盖A, 这里Q=P1...P_{n-2}
%是一系列%Householder矩阵的积且Q^T(A-aI)(A-bI)是上三角形矩阵.
I3=eye(3); I2=eye(2); s=A(n-1,n-1)+A(n,n);
t=A(n-1,n-1)*A(n,n)-A(n-1,n)*A(n,n-1);
x=A(1,1)^2+A(1,2)*A(2,1)-s*A(1,1)+t;
y=A(2,1)*(A(1,1)+A(2,2)-s); z=A(2,1)*A(3,2);
for k=0:n-3
   [v,beta]=mhouse([x,y,z]');  q=max(1,k); 
   A(k+1:k+3,q:n)=(I3-beta*v*v')*A(k+1:k+3,q:n);
   r=min(k+4,n);
   A(1:r,k+1:k+3)=A(1:r,k+1:k+3)*(I3-beta*v*v');
   x=A(k+2,k+1); y=A(k+3,k+1);
   if (k<n-3), z=A(k+4,k+1); end
end
[v,beta]=mhouse([x,y]');
A(n-1:n,n-2:n)=(I2-beta*v*v')*A(n-1:n,n-2:n);
A(1:n,n-1:n)=A(1:n,n-1:n)*(I2-beta*v*v');
