%ʵ������Householder�任����--r_house.m
function [v,beta]=r_house(x)
%����������Householder����H=I-beta*v*v'������v(1)=1��v��beta.
n=length(x);
eta=norm(x,inf); x=x/eta;
sigma=x(2:n)'*x(2:n);
v=[1; x(2:n)];
if sigma==0
     if x(1)>=0
         beta=0;
     else
         beta=2;
     end
else
    alpha=(x(1)^2+sigma)^0.5;
    if x(1)<=0
        v(1)=x(1)-alpha;
    else
        v(1)=-sigma/(x(1)+alpha);
    end
    beta=2*v(1)^2/(sigma+v(1)^2);
    v=v/v(1);
end
 