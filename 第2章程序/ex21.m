%Àı2.1-ex21
x=[2,5,7,1]'; [v,beta]=r_house(x);
H=eye(length(x))-beta*v*v';
y=H*x