%ex438
clear all
n=15; e=[1:n]'; e1=ones(n-1,1); e2=ones(n-2,1); 
e3=ones(n-3,1); e4=ones(n-4,1);
A=diag(e)+diag(e1,1)-diag(e1,-1)-2*diag(e2,2)+2*diag(e2,-2);
A=A+3*diag(e3,3)-3*diag(e3,-3)-2*diag(e4,4)+2*diag(e4,-4);
[Lam, V, iter, ki]=ddiqrm_vec(A);
D=eig(A); %��ϵͳ������A��ȫ������ֵ
disp([Lam,  D])  %��ʾ���
iter   %��ʾ��������
ki   %��ʾ���ݷ���ÿ�����������ĵ�������
err=norm(A*V-V*diag(Lam),inf)  %��֤�����Ƿ���ȷ
