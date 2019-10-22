%��������Householder�任����--house.m
function [w,gamma]=c_house(x)
%����������x, ����������һ������w��һ����gamma
%����H*x=gamma*e1,����H=I-w*w',||w||=sqrt(2).
w=x; gamma=norm(x);  
if gamma==0
    w(1)=sqrt(2); return;
end
if w(1)==0
    tau=1;
else
    tau=conj(w(1))/abs(w(1));
end
w=(tau/gamma)*w; w(1)=w(1)+1;
w=w/sqrt(w(1)); 
gamma=-conj(tau)*gamma;
    
 