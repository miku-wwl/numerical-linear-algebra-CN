function [M,Q,iter] =rest_co_Lanczos(A,m)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
m=m+2; k=2*m; Q=[ ]; tol=1.0e-5;
v=rand(size(A,1),1);
q=v/norm(v); Q=[Q,q];
u=A*q; alpha(1)=u'*q;
u=u-alpha(1)*q; u=u-(u'*q)*q;
beta(1)=norm(u); q=u/beta(1);
Q=[Q,q]; 
[Q,alpha,beta] = expand_lanczos(A,Q,alpha,beta,1,m);
iter=0;
while(iter<100)
    iter=iter+1;
    [Q,alpha,beta] = expand_lanczos(A,Q,alpha,beta,m,k);
    T=diag(alpha)+diag(beta(1:end-1),1)+diag(beta(1:end-1),-1);
    [Y,M]=eig(T);  %M的对角元按升幂按列
    if beta(end)*abs(Y(k,k))<tol
       break; 
    end
    [M,I]=sort(-diag(M));
    M=-diag(M);  Y=Y(:,I);
    Q=Q*Y(:,1:m); z=Y(k,1:m); M=M(1:m,1:m);
    [V,alpha,beta,theta] = trireduce(M,z);
    Q=Q*V'; beta(m)=theta*beta(k); Q(:,m+1)=Q(:,k+1);
    T=V*M*V';
end




