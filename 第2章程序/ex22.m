%Àý2.2-ex22
x=[2+i,5-3*i,7,1+2*i]'; 
[w,gama]=c_house(x);
H=eye(length(x))-w*w'; 
y=H*x
