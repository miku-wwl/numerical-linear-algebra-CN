function [Iter,time]=ex610(alpha) 
%clc;  
A=[2 2 1; 2 5 4;1 4 5]; B=[4 1 0; 1 4 1;0 1 3];
C=[2 1 0; 1 2 2; 0 2 3]; D=[2 0 1; 0 2 2; 1 2 3];
F=eye(3);
A1=alpha*A+C; C1=alpha*A-C;
D1=alpha*D+B; B1=alpha*D-B;
rA1=inv(A1);   rD1=inv(D1);
P=rA1*C1;   Q=B1*rD1;
Y0=2*alpha*rA1*F*rD1';
X0=Y0; Iter=0;
tic
while(Iter<=500)
    X=P*X0*Q+Y0;
    Err=norm(X-X0,'fro');
    if(Err<=1.e-10)
       break;
    end
    X0=X;  Iter=Iter+1;
end
time=toc;
