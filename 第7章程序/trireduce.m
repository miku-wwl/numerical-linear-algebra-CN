function [Q,alpha,beta,tau]=trireduce(M,z)
%给定m阶对角矩阵M和m维向量z,本算法计算正交矩阵Q和对
%称三对角矩阵T及实数alpha,使得Q*z=tau*e_m,QMQ'=T
m=size(M,1);
for i=1:m
    alpha(i)=M(i,i);
end
beta=zeros(m-1,1); Q=eye(m);
for i=1:m-1
    delta=sqrt(z(i)^2+z(i+1)^2); 
    c=z(i+1)/delta; s=-z(i)/delta;
    z(i+1)=delta; z(i)=0;
    a=[c,s;-s,c]*[alpha(i),beta(i); beta(i), alpha(i+1)]*[c,s;-s,c]';
    alpha(i)=a(1,1); beta(i)=a(1,2); alpha(i+1)=a(2,2); 
    Q(i:i+1,:)=[c,s;-s,c]*Q(i:i+1,:);
    if i>1
        gamma=-s*beta(i-1);  beta(i-1)=c*beta(i-1);      
        for r=i:-1:2
            sigma=sqrt(gamma^2+beta(r)^2);
            c=beta(r)/sigma; s=-gamma/sigma;
            a=[c,s;-s,c]*[alpha(r-1),beta(r-1); beta(r-1), alpha(r)]*[c,s;-s,c]';
            alpha(r-1)=a(1,1); beta(r-1)=a(1,2); alpha(r)=a(2,2);
            beta(r)=sigma;
            if r>2
                 gamma=-s*beta(r-2);  beta(r-2)=c*beta(r-2);   
            end
            Q(r-1:r,:)=[c,s;-s,c]*Q(r-1:r,:);
        end
    end
end
tau=z(m);  %T=M;
% T=diag(alpha)+diag(beta,1)+diag(beta,-1);
% norm(T-V*M*V')
            
            


