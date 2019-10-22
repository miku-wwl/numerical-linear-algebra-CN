function [M,V,iter] =irest_Lanczos(A,m)
%给定一个大型稀疏的对称矩阵A,本算法计算A的m个两端特征对(m<<n).
m=m+2; k=2*m; V=[ ];  tol=1.0e-30; iter=0;
v=rand(size(A,1),1); v=v/norm(v); V=[V,v];
u=A*v; alpha(1)=u'*v;
u=u-alpha(1)*v; u=u-(u'*v)*v;
beta(1)=norm(u); v=u/beta(1); V=[V,v]; 
[V,alpha,beta]=expand_lanczos(A,V,alpha,beta,1,m);
while(iter<1000)
    iter=iter+1;
    [V,alpha,beta]=expand_lanczos(A,V,alpha,beta,m,k);
    betak=beta(end);  %将beta的最后一个分量存起来
    vk1=V(:,end);   %将V的最后一列存起来
    T=diag(alpha)+diag(beta(1:end-1),1)+diag(beta(1:end-1),-1);
    [Y,M]=eig(T);  %M的对角元按升序按列
    [M,I]=sort(diag(M),'descend'); %对M的对角元按降序按列
    M=diag(M);  Y=Y(:,I); %第k列到了第1列
    V=V(:,1:k)*Y(:,1:m); z=Y(k,1:m)'; M=M(1:m,1:m);
    if abs(beta(end))*abs(Y(k,1))<tol
       break; 
    end
    [Q,alpha,beta,theta] = trireduce2(M,z);
    V=V*Q';  V=[V,vk1];  
    beta(m)=theta*betak;   
end





