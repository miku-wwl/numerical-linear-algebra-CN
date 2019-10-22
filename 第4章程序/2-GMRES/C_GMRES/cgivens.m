%Givens rotation matrix 
function [c, s, eta] = cgivens(a,b)
if b==0
    c=1;s=0; eta=a;
end
if a==0
    c=0; s=1; eta=b;
end
mu=a/abs(a); tau=abs(a)+abs(b);
delta=tau*sqrt(abs(a/tau)^2+abs(b/tau)^2);
c=abs(a)/delta; s=mu*conj(b)/delta; eta=delta*mu;
        
        