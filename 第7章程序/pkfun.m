function p=pkfun(t,a,b)
%本函数计算Sturm序列{p_k(t)}的值
n=length(a);
p(1)=a(1)-t; p(2)=(a(2)-t)*p(1)-b(1);
for k=3:n
    p(k)=(a(k)-t)*p(k-1)-b(k-1)^2*p(k-2);
end
