function q=qkfun(t,a,b)
%�����������ֹ����Ķ���ʽ����{q_k(t)}��ֵ
n=length(a);   
q(1)=a(1)-t; 
for k=2:n
    if q(k-1)==0,
        q(k)=1;
    else
        q(k)=a(k)-t-b(k-1)^2/q(k-1);
    end
end
