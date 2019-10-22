%¸´Given±ä»»³ÌÐò: c_givens.m
function [c,s,eta]=c_givens(a,b)
%¼ÆËã c,s, Âú×ã[c s; -s' c]*[a;b]=[eta;0]
if b==0,  c=1; s=0; eta=a; end
if a==0,  c=0; s=1; eta=b; end
u=a/abs(a); t=sqrt(abs(a)^2+abs(b)^2);
c=abs(a)/t; s=u*conj(b)/t; eta=u*t;

