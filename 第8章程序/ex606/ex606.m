function [Iter, time]=ex606(alpha) 
%clc;  %调用时要注意alpha是负数
ni=3; m=100;  I=eye(ni*m); F=I;
A = cell(m,m);  AI = [4 -1 -1; -1 4 -1; 0 -1  4];  Ii = -eye(ni);
tic
for i=1:m
    for j=1:m
        A{i,j}=zeros(ni);
    end
end
for i=2:m-1
     A{i,i-1}=Ii; A{i,i}=AI; A{i,i+1}=Ii;
end
A{1,1}=AI;      A{1,2}=Ii; 
A{m,m-1}=Ii; A{m,m}=AI;
A=cell2mat(A);
%l1=max(eig(A))  %最大特征值
%ln=min(eig(A))   %最小特征值
%aopt=-1/sqrt(l1*ln)     %最优参数
B=I+2*inv(alpha*A-I);
Y0=-0.5*alpha*(B-I)*F*(B-I)';
X0=Y0; Iter=0;
while(Iter<=500)
    X=B*X0*B'+Y0;
    Err=norm(X-X0,'fro');
    if(Err<=1.e-10)
       break;
    end
    X0=X;  Iter=Iter+1;
end
time=toc;
