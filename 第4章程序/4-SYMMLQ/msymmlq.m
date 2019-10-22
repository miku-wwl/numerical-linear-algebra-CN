%SYMMLQ方法程序-msymmlq.m
function [x,k,time,res,resvec]=msymmlq(A,b,x,max_it,tol)
%SYMMLQ求解对称不定方程组Ax=b
tic; r=b-A*x;  beta1=norm(r); v1=r/beta1;
z=A*v1; alpha=v1'*z; v_hat=z-alpha*v1; 
beta=norm(v_hat); v2=v_hat/beta; w_hat=v1;
nu_hat=alpha; delta_hat=beta;
zeta=0; eta=0; 
nu1=sqrt(alpha^2+beta^2);
c1=nu_hat/nu1; s1=beta/nu1;
w1=c1*w_hat+s1*v2; w_hat=s1*w_hat-c1*v2;
zeta1=beta1/nu1;  tx=x+zeta1*w1;
k=1; res=1;
while (k<max_it)   
    res=res*s1;   resvec(k)=res;
    if (res<tol),   break;  end
    z=A*v2; alpha=v2'*z;  
    v_hat=z-alpha*v2-beta*v1; 
    beta=norm(v_hat); v3=v_hat/beta;
    nu_hat=s1*delta_hat-c1*alpha;
    nu2=sqrt(nu_hat^2+beta^2);
    c2=nu_hat/nu2; s2=beta/nu2;
    delta2=c1*delta_hat+s1*alpha;
    zeta2=-(eta*zeta+delta2*zeta1)/nu2;
    zeta_hat=zeta2/c2;
    w2=c2*w_hat+s2*v3; 
    x=tx+zeta_hat*w_hat;  tx=tx+zeta2*w2;
    delta_hat=-c1*beta; eta=s1*beta;
    w_hat=s2*w_hat-c2*v3;
    k=k+1; v1=v2; v2=v3; c1=c2; s1=s2; 
    zeta=zeta1; zeta1=zeta2;
end
time=toc;

