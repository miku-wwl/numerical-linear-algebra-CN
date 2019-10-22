a=3+4*i; b=2+5*i;
[c,s,eta]=c_givens(a,b);
[c,s;-s', c]*[a, b]'
eta