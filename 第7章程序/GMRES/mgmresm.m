%重开始GMRES(m)
function [x, out, int, time, res, resvec, flag] = mgmresm( A, b, x, restrt, max_it, tol )
flag = 0; int=0; bnrm2 = norm(b);
if  (bnrm2 == 0.0), bnrm2 = 1.0; end
r =  b-A*x;  %计算残差
res= norm(r) / bnrm2; resvec(1)=res;
if (res < tol ) return, end
n =length(b);  m = restrt;
Q(1:n,1:m+1) = zeros(n,m+1);
H(1:m+1,1:m) = zeros(m+1,m);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1 = zeros(n,1);
e1(1) = 1.0;
%
tic;
for k = 1:max_it,      % 迭代开始
    Q(:,1) = r / norm( r );
    xi = norm( r )*e1;
    for j = 1 : m,          %用Arnoldi方法构造正交基 
        w = A*Q(:,j);     
        for i = 1 : j
            H(i,j)= w'*Q(:,i);
            w = w - H(i,j)*Q(:,i);
        end
        H(j+1,j) = norm(w);
        if abs(H(j+1,j))/bnrm2<tol,
            return;
        else
            Q(:,j+1) = w / H(j+1,j);
        end
        for i = 1:j-1,                              % apply Givens rotation
            temp     =  cs(i)*H(i,j) + sn(i)*H(i+1,j);
            H(i+1,j) = -sn(i)*H(i,j) + cs(i)*H(i+1,j);
            H(i,j)   = temp;
        end
        [cs(j),sn(j),H(j,j)] = givens( H(j,j), H(j+1,j) ); % form j-th Givens rotation matrix
        % approximate residual norm
        xi(j+1) = -sn(j)*xi(j);
        xi(j)   = cs(j)*xi(j);
        %H(j,j) = cs(j)*H(j,j) + sn(j)*H(j+1,j);
        H(j+1,j) = 0.0;
        res  = abs(xi(j+1)) / bnrm2;
        resvec((k-1)*m+j+1)=res;
        if (res <= tol ),                         
            y = H(1:j,1:j) \xi(1:j);               
            x = x + Q(:,1:j)*y; 
            break; %跳出内循环
        end
    end
    if (res <= tol )
        out=k; int=j; break;   %跳出外循环
    end
    y = H(1:m,1:m) \xi(1:m);
    x = x + Q(:,1:m)*y;                          
    r = b-A*x ;       
end
if (res > tol ) flag = 1; end;          % not converged
time=toc;
% END of gmres.m