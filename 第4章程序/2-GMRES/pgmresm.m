%�ؿ�ʼPGMRES����-PGMRES(m)-gmresm.m
function [x,out,int,time,res,resvec,flag]=pgmresm(A,b,x,M,restrt,max_it,tol )
tic; flag=0; int=0; r=M\(b-A*x);  %����в�
beta=norm(r); res=norm(r)/beta; resvec(1)=res;
n =length(b);  m=restrt;
V(1:n,1:m+1)=zeros(n,m+1);
H(1:m+1,1:m)=zeros(m+1,m);
c(1:m)=zeros(m,1); s(1:m)=zeros(m,1);
e1=zeros(n,1); e1(1)=1.0;
for k=1:max_it,  
    V(:,1)=r/norm(r); xi=norm(r)*e1;
    for j=1:m,  %��Arnoldi�������������� 
        w=M\(A*V(:,j));     
        for i =1:j
            H(i,j)=w'*V(:,i);  w=w-H(i,j)*V(:,i);
        end
        H(j+1,j)=norm(w);
        if abs(H(j+1,j))/beta<tol,
            return;
        else
            V(:,j+1)=w/H(j+1,j);
        end
        for i=1:j-1, %��i��Givens�任
            temp=c(i)*H(i,j)+s(i)*H(i+1,j);
            H(i+1,j)=-s(i)*H(i,j)+c(i)*H(i+1,j);
            H(i,j)=temp;
        end
        [c(j),s(j),H(j,j)]=givens(H(j,j),H(j+1,j) ); %��j��Givens�任
        xi(j+1)=-s(j)*xi(j); xi(j)=c(j)*xi(j); H(j+1,j)=0.0;
        res=abs(xi(j+1))/beta; resvec((k-1)*m+j+1)=res;
        if (res<=tol),                         
            y=H(1:j,1:j)\xi(1:j); x=x+V(:,1:j)*y; 
            break; %������ѭ��
        end
    end
    if (res<tol )
        out=k; int=j; break;  %������ѭ��
    end
    y=H(1:m,1:m)\xi(1:m);
    x=x+V(:,1:m)*y;                          
    r=M\(b-A*x);       
end
if (res>tol), flag=1; end;  %������
time=toc;