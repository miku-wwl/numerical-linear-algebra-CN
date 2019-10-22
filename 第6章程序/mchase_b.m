%��׷�Ϸ�����--mchase_b.m
function [x]=mchase_b(Ai,Bi,Ci,fi,m)
%��׷�Ϸ������ԽǷ����� Ax=f
%����: AiΪA�Ĵ��¶Խǿ�,BiΪA�����Խǿ�,
%CiΪA�Ĵ��϶Խǿ�, fiΪ�Ҷ�����. 
%���: ������x (LU�ֽ��е�l(k),u(k)�����b(k),
%c(k)��λ��, y(k)��x(k)�Ⱥ�����d(k)��λ��)
x=cell(m,1);
L{1}=Bi{1}; y{1}=L{1}\fi{1};
for k=2:m  
    U{k-1}=L{k-1}\Ci{k-1};
    L{k}=Bi{k}-Ai{k}*U{k-1};
    y{k}=L{k}\(fi{k}-Ai{k}*y{k-1});
end  
x{m}=y{m};
for k=m-1:-1:1
    x{k}=y{k}-U{k}*x{k+1};
end
    