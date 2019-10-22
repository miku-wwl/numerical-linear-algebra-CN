function [c,s,gama]=givens(a,b)
%实数情形的Givens变换
%计算 c=cos(theta),s=sin(theta), 满足[c -s; s c]*[a;b]=[gama;0]
if b==0
    c=1;s=0; gama=a;
end
if a==0
    c=0; s=1; gama=b;
end
if abs(b)>=abs(a)
    t=a/b; s=1/sqrt(1+t^2); c=s*t;
else
    t=b/a;  c=1/sqrt(1+t^2); s=c*t;
end
gama=c*a+s*b;

