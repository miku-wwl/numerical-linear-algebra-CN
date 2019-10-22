function [M,V,iter] =irest_Lanczos(A,m)
%����һ������ϡ��ĶԳƾ���A,���㷨����A��m������������(m<<n).
m=m+2; k=2*m; V=[ ];  tol=1.0e-30; iter=0;
v=rand(size(A,1),1); v=v/norm(v); V=[V,v];
u=A*v; alpha(1)=u'*v;
u=u-alpha(1)*v; u=u-(u'*v)*v;
beta(1)=norm(u); v=u/beta(1); V=[V,v]; 
[V,alpha,beta]=expand_lanczos(A,V,alpha,beta,1,m);
while(iter<1000)
    iter=iter+1;
    [V,alpha,beta]=expand_lanczos(A,V,alpha,beta,m,k);
    betak=beta(end);  %��beta�����һ������������
    vk1=V(:,end);   %��V�����һ�д�����
    T=diag(alpha)+diag(beta(1:end-1),1)+diag(beta(1:end-1),-1);
    [Y,M]=eig(T);  %M�ĶԽ�Ԫ��������
    [M,I]=sort(diag(M),'descend'); %��M�ĶԽ�Ԫ��������
    M=diag(M);  Y=Y(:,I); %��k�е��˵�1��
    V=V(:,1:k)*Y(:,1:m); z=Y(k,1:m)'; M=M(1:m,1:m);
    if abs(beta(end))*abs(Y(k,1))<tol
       break; 
    end
    [Q,alpha,beta,theta] = trireduce2(M,z);
    V=V*Q';  V=[V,vk1];  
    beta(m)=theta*betak;   
end





