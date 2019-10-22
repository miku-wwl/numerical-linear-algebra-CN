%����Givens-Householder�����������������ķ�����������-ghvector.m
function  [Lam,V,ki]=ghvector(A,tol,max_it)
%�������÷���������ʵ�Գƾ���A�ĵ�m�����������.
%�������:AΪn�׶ԳƷ���,mΪ��������������ֵ���,tolΪ�������.
%�������:lambdaΪ��m�������ֵ,xΪ��Ӧ����������,kΪ��������.
n=size(A,1);  Lam=zeros(n,1); x=rand(n,1); %x=ones(n,1);
[T,Q]=mhessen(A); %������Hessenberg������
for i=1:n
    [la]=givens_househ(T,i,tol,max_it); Lam(i)=la;
    [lam,v,k2]=mvpower(T,x,Lam(i)); %���÷��ݷ�����
    Lam(i)=lam; V(:,i)=v; ki(i)=k2; 
    %norm(T*v-Lam(i)*v),
end
V=Q*V; %V��ÿһ��Ϊ��������

    

   