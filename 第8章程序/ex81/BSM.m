%Bartels-Rtewart Method- AX+XB=F
function [X]=BSM(A,B,F)
[m,n]=size(F); X=zeros(m,n);
[U,R]=schur(A);
if (norm(B-A,'fro')==0)
    V=U; S=R;
else if (norm(B-A','fro')==0)
        id=[size(A,1):-1:1];
        V=U(:,id); S=R(id,id)';
    else
        [V,S]=schur(B);
    end
end
F=U'*F*V; Rsq=R*R; I=eye(m); k=1;
while (k<n+1)
    if (k<n-1)&(abs(S(k+1,k))>10*eps*max(abs(S(k,k)),abs(S(k+1,k+1))))
        s11=S(k,k);     s12=S(k,k+1);
        s21=S(k+1,k);  s22=S(k+1,k+1);
        b=F(:,k:k+1)-X(:,1:k-1)*S(1:k-1,k:k+1);
        b=[R*b(:,1)+s22*b(:,1)-s21*b(:,2), R*b(:,2)+s11*b(:,2)-s12*b(:,1)];
        x=(Rsq+(s11+s22)*R+(s11*s22-s21*s12)*I)\b;
        X(:,k:k+1)=x;   k=k+2;
    else
        b=F(:,k)-X(:,1:k-1)*S(1:k-1,k);
        x=(R+S(k,k)*I)\b; 
        X(:,k)=x;   k=k+1;
    end
end
X = U*X*V';


    

