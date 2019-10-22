%GMRES法
function [x, k, time, res, resvec, flag] = mgmres( A, b, x0, max_it, tol )
flag = 0;  bnrm2 = norm(b);
if  (bnrm2 == 0.0), bnrm2 = 1.0; end
r =  b-A*x0;  %计算残差
res= norm(r) / bnrm2; resvec(1)=res;
if (res < tol ) return, end
n =length(b);  e1 = zeros(n,1); e1(1) = 1.0;
%
tic;
Q(:,1) = r / norm( r );
xi = norm( r )*e1; k=0;
while (k <= max_it)   
    k=k+1;
    w = A*Q(:,k);   
    for i = 1 : k
        H(i,k)= w'*Q(:,i);
        w = w - H(i,k)*Q(:,i);
    end
    H(k+1,k) = norm(w);
    if abs(H(k+1,k))/bnrm2<tol,
        return;
    else
        Q(:,k+1) = w / H(k+1,k);
    end
    for i = 1:k-1,                            
        temp     =  c(i)*H(i,k) + s(i)*H(i+1,k);
        H(i+1,k) = -s(i)*H(i,k) +c(i)*H(i+1,k);
        H(i,k)   = temp;
    end
    [c(k),s(k),H(k,k)] = givens(H(k,k), H(k+1,k) ); % form k-th Givens rotation matrix
    % approximate residual norm
    xi(k+1) = -s(k)*xi(k);
    xi(k)  = c(k)*xi(k);     
    %H(k,k) = c(k)*H(k,k) + s(k)*H(k+1,k);
    H(k+1,k) = 0.0;
    res  = abs(xi(k+1)) / bnrm2;
    resvec(k+1)=res;
    if (res <= tol ),   
        y = H(1:k,1:k) \xi(1:k); 
        x = x0 + Q(:,1:k)*y; 
        break; %跳出循环
    end
end
y = H(1:k,1:k) \xi(1:k);    x = x0 + Q(:,1:k)*y; 
if (res > tol ) flag = 1; end;          % not converged
time=toc;
% END of gmres.m