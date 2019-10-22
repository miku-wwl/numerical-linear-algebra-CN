%基于Givens-Householder方法的求特征向量的反迭代法程序-ghvector.m
function  [Lam,V,ki]=ghvector(A,tol,max_it)
%本程序用反迭代法求实对称矩阵A的第m大的特征向量.
%输入参数:A为n阶对称方阵,m为按降幂排列特征值序号,tol为容许误差.
%输出参数:lambda为第m大的特征值,x为对应的特征向量,k为迭代次数.
n=size(A,1);  Lam=zeros(n,1); x=rand(n,1); %x=ones(n,1);
[T,Q]=mhessen(A); %调用上Hessenberg化程序
for i=1:n
    [la]=givens_househ(T,i,tol,max_it); Lam(i)=la;
    [lam,v,k2]=mvpower(T,x,Lam(i)); %调用反幂法程序
    Lam(i)=lam; V(:,i)=v; ki(i)=k2; 
    %norm(T*v-Lam(i)*v),
end
V=Q*V; %V的每一列为特征向量

    

   