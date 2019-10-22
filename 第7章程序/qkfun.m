function q=qkfun(t,a,b)
%本函数计算防止溢出的多项式序列{q_k(t)}的值
n=length(a);   
q(1)=a(1)-t; 
for k=2:n
    if q(k-1)==0,
        q(k)=1;
    else
        q(k)=a(k)-t-b(k-1)^2/q(k-1);
    end
end
