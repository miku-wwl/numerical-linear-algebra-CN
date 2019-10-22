function [v,beta]=mhouse(x)
%����nά����x,��������������v(1)=1������v����beta,
%����betaʹ��P=I-beta*v*v'������������Px=||x||e_1.
n=length(x); s=x(2:n)'*x(2:n); v=[1;x(2:n)];
if s<=eps
    beta=0;
else
    mu=sqrt(x(1)^2+s);
    if x(1)<=0
        v(1)=x(1)-mu;
    else
        v(1)=-s/(x(1)+mu);
    end
    beta=2*v(1)^2/(s+v(1)^2);
    v=v/v(1);
end