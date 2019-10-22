function [x, out, int, time, res, resvec, flag] = cpgmresm( A, b, x, M, restrt, max_it, tol )
flag = 0; bnrm2 = norm(M\b);
if  (bnrm2 == 0.0), bnrm2 = 1.0; end
r = M \ ( b-A*x );  %计算残差r=M^{-1}r0
res = norm(r) / bnrm2; resvec(1)=res;
if (res < tol), return, end
n= length(b);   m = restrt;
Q(1:n,1:m+1) = zeros(n,m+1);
H(1:m+1,1:m) = zeros(m+1,m);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1 = zeros(n,1); e1(1) = 1.0;
%
tic;
for k = 1:max_it,      % 迭代开始
    r = M \ ( b-A*x );
    Q(:,1) = r / norm( r );
    xi = norm( r )*e1;
    for j = 1 : m,              %用Arnoldi方法构造正交基 
        w = M \ (A*Q(:,j));     
        for i = 1 : j
            H(i,j)= conj(w')*conj(Q(:,i));
            w = w - H(i,j)*Q(:,i);
        end
        H(j+1,j) = norm(w);
        Q(:,j+1) = w / H(j+1,j);
        for i = 1:j-1,               % apply Givens rotation
            temp     =  cs(i)*H(i,j) + sn(i)*H(i+1,j);
            H(i+1,j) = -conj(sn(i))*H(i,j) + conj(cs(i))*H(i+1,j);
            H(i,j)   = temp;
        end
        [cs(j),sn(j),H(j,j)] = cgivens( H(j,j), H(j+1,j) ); % form j-th Givens rotation matrix
        % approximate residual norm
        xi(j+1) = -conj(sn(j))*xi(j);
        xi(j)   =  cs(j)*xi(j);  
        %H(j,j) = cs(j)*H(j,j) + sn(j)*H(j+1,j);
        H(j+1,j) = 0.0;
        res  = abs(xi(j+1)) / bnrm2;
        resvec((k-1)*m+j+1)=res;
        if (res <= tol ),                         
            y = H(1:j,1:j) \xi(1:j);               
            x = x + Q(:,1:j)*y;
            break;
        end
    end
    if (res <= tol )  
        out=k; int=j; break;
    end
    y = H(1:m,1:m) \xi(1:m);
    x = x + Q(:,1:m)*y;                          
    r = M \ ( b-A*x );       
end
if (res > tol ) flag = 1; end;    % not converged
time=toc;
% END of gmres.m