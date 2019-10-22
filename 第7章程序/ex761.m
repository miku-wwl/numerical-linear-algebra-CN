%例7.6.1-ex761.m
clear all
n=500; v=rand(n,1); k=40;
e=randn(n,1); A=diag(e);
ev=eig(A);
[mu,U,lambda] = Class_Lanczos(A,v,k);
plot(size(mu,2)+1,ev,'+k');  hold on
for i=1:size(mu,2)
    plot(i,mu{i},'+k');
end
hold off
xlabel('子空间维数k'), ylabel('近似特征值'),
%axis([0,k+2,-3,3])





