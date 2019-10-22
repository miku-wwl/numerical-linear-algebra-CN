%QMRES·½·¨³ÌÐò-qmres.m
function [x,k,time,res,resvec]=mqmres(A,b,x,max_it, tol)
tic; r=b-A*x;  bt=norm(r); v1=r/bt; 
w1=v1;  av=A*v1; atw=A'*w1; alpha=w1'*av;
v=av-alpha*v1;  w=atw-alpha*w1; omega=w'*v;
if (abs(omega)<tol)
    x=x+r/alpha; return;
else
    beta=sqrt(abs(omega)); gama=omega/beta;
    v2=v/beta; w2=w/gama;
end
[c1,s1,sigma]=givens(alpha,beta);
p1=v1/sigma;  rho=-bt*s1; tau=bt*c1;
p0=zeros(length(b),1); c0=1; s0=0; 
x=x+tau*p1;  k=1;  resvec(1)=1;
while (k<max_it)
    av=A*v2; atw=A'*w2;  alpha=w2'*av; 
    v=av-alpha*v2-gama*v1;  w=atw-alpha*w2-beta*w1;  
    omega=w'*v;
    if (abs(omega)==0)
        break;
    else
        beta=sqrt(abs(omega)); gama1=omega/beta;
        v3=v/beta; w3=w/gama1;
    end
    epsi=s0*gama; hgama=c0*gama;
    delta=c1*hgama+s1*alpha; halpha=-s1*hgama+c1*alpha;
    [c2,s2,sigma]=givens(halpha,beta);
    tau=rho*c2; rho=-rho*s2;
    p2=(v2-epsi*p0-delta*p1)/sigma;
    x=x+tau*p2; 
    res=abs(rho)/bt; resvec(k+1)=res;
    if (res<tol), break; end
    v1=v2; v2=v3; w1=w2; w2=w3; 
    p0=p1; p1=p2; gama=gama1; 
    c0=c1; c1=c2; s0=s1; s1=s2;
    k=k+1;
end
time=toc;

