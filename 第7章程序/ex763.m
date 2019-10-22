%Àý7.6.3
clear all
n=50; v=ones(n,1); k=30;
P=rand(n,n); A1=diag(-24:25);
A=(P\A1)*P;
[mu,U] = Class_Arnoldi(A,v,k);
ev=eig(A);
for i=1:k
    plot(i,mu{i},'+k'); hold on
end
plot(k+1,ev,'+k')
hold off
axis([0,33,-25,25,])





