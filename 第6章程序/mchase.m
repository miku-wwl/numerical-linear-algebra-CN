%׷�Ϸ�����--mchase.m
function [f]=mchase(a,b,c,f)
%��׷�Ϸ������ԽǷ����� Ax=f
%����: aΪA���¶Խ���,bΪA�����Խ�
%�� cΪA�Ĵ��϶Խ���, fΪ�Ҷ�����. 
%���: ������f (LU�ֽ��е�l(k),u(k)�����b(k),
%c(k)��λ��, y(k)��x(k)�Ⱥ�����d(k)��λ��)
n=length(b); f(1)=f(1)/b(1);
for k=2:n  
    c(k-1)=c(k-1)/b(k-1);
    b(k)=b(k)-a(k)*c(k-1);
    f(k)=(f(k)-a(k)*f(k-1))/b(k);
end    
for k=n-1:-1:1
    f(k)=f(k)-c(k)*f(k+1);
end
    