function [rB, Iter, Err]=ex605(alpha) 
%clc;  %调用时要注意alpha是负数
A = [4 -1 -1; -1 4 -1; -1 -1  4]; 
I = eye(3);  F = eye(3);
B = I + 2*inv(alpha*A-I);
rB=max(abs(eig(B)));
Y0=-0.5*alpha*(B-I)*F*(B-I)';
X0=Y0; Iter=0;
while(Iter<=500)
    X=B*X0*B'+Y0;
    Err=norm(X-X0,'fro');
    if(Err<=1.e-6)
       break;
    end
    X0=X;  Iter=Iter+1;
end
