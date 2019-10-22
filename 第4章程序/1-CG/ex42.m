%��4.2
clear all
n=1000; tol=1e-10;  max_it=1000;
e=[1:n]'; e1=ones(n-1,1);
A=diag(e)-diag(e1,1)-diag(e1,-1);
xs=ones(n,1);  b=A*xs;  x0=zeros(n,1);
M=diag(diag(A));
[x1,iter1,time1, res1, resvec1] = mcg(A,b,x0,max_it,tol);
[x2,iter2,time2, res2, resvec2] = pcg(A,b,x0,M,max_it,tol);
%��ʽ����ʾ
fid = 1;
fprintf(fid, ' �� ��     Iter         CPU         RES      \n' );     fprintf('\n');
fprintf(fid, 'PCG      %3i       %8.4f    %11.4e\n', iter1, time1, res1); fprintf('\n');
fprintf(fid, ' CG       %3i       %8.4f    %11.4e\n', iter2, time2, res2); fprintf('\n');
%���ӻ�
t=1:length(resvec1); t2=1:length(resvec2);
semilogy(t,resvec1,'.-k'); hold on
semilogy(t2,resvec2,'*-k'); hold off
legend('CG ����','PCG ����')
er1=norm(x1-xs), err2=norm(x2-xs),
rs1=norm(b-A*x1), rs2=norm(b-A*x2)
xlabel('��������'); ylabel('��Բв��');
%axis([1,iter2,1.e-12, 1.e0])

