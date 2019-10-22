%׷�Ϸ�����--mchase.m
function [x]=mchaseold(a,b,c,f)
%��׷�Ϸ������ԽǷ����� Ax=f
%����: a=[0, a(2),a(3),...,a(n)] ΪA���¶Խ��ߣ�
%b=[b(1),b(2),...,b(n-1),b(n)] ΪA�����Խ��ߣ�
%c=[c(1),c(2),...,c(n-1),0] ΪA�Ĵ��϶Խ���, 
%fΪ�Ҷ�����. ���: ������x
n=length(b);
l(1)=b(1); y(1)=f(1)/l(1);
for k=2:n  
    u(k-1)=c(k-1)/l(k-1);
    l(k)=b(k)-a(k)*u(k-1);
    y(k)=(f(k)-a(k)*y(k-1))/l(k);
end    
x(n)=y(n);
for k=n-1:-1:1
    x(k)=y(k)-u(k)*x(k+1);
end
    