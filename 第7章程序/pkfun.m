function p=pkfun(t,a,b)
%����������Sturm����{p_k(t)}��ֵ
n=length(a);
p(1)=a(1)-t; p(2)=(a(2)-t)*p(1)-b(1);
for k=3:n
    p(k)=(a(k)-t)*p(k-1)-b(k-1)^2*p(k-2);
end
