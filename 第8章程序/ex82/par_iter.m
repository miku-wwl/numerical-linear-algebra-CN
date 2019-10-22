function [Iter, time]=par_iter(A,B,F,alpha,tol,max_it) 
if nargin<6, max_it=1000;end
if nargin<5,tol=1e-10;end
if nargin<4,alpha=-1;end
tic, I=eye(size(A));
C=I+2*((alpha*A-I)\I);
D=I+2*((alpha*B-I)\I);
Y0=-0.5*alpha*(C-I)*F*(D-I);
X0=Y0;  Iter=0;
while(Iter<=max_it)
    X=C*X0*D+Y0;
    Err=norm(X-X0,'fro');
    if(Err<=tol), break; end
    X0=X;  Iter=Iter+1;
end
time=toc;
