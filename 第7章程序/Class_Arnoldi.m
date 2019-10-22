function [mu,U] = Class_Arnoldi(A,v,k)
%æ≠µ‰Arnoldi∑Ω∑®
tol=1.e-12;
q1=v/norm(v); Q(:,1)=q1;
for k1=1:k
    [Q,H] = Arnoldi2(A,q1,k1);
    [Y,Mu]=eig(H(1:k1,1:k1));
    mu{k1}=diag(real(Mu));
    U=Q(:,1:k1)*Y;
    if abs(H(k1+1,k1))*abs(Y(k1,k1))<tol
        break;
    end
end



