%Àý6.14-ex614.m
clear all; n=128;  e=ones(n-1,1);
A=2*eye(n)+diag(e,-1)+diag(e,1);
tic, [kappa]=cond_inf(A), toc
tic, conds=norm(A,inf)*norm(inv(A),inf), toc