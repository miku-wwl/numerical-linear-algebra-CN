%˫����������--mdouble_par.m
function [x]=mdouble_par(Ai,Bi,Ci,fi,m)
%��˫�����������ԽǷ����� Ax=f
%����: AiΪA�Ĵ��¶Խǿ�,BiΪA�����Խǿ�,
%CiΪA�Ĵ��϶Խǿ�, fiΪ�Ҷ�����.
%���: ������x
x=cell(m,1); s=cell(m,1); T=cell(m,1);
s{1}=zeros(3,1); s{2}=Ci{1}\fi{1}; 
T{1}=eye(3); T{2}=-Ci{1}\Bi{1};
for k=2:m-1 
    s{k+1}=Ci{k}\(fi{k}-Ai{k}*s{k-1}-Bi{k}*s{k});
    T{k+1}=-Ci{k}\(Ai{k}*T{k-1}+Bi{k}*T{k});
end  
x{1}=(Ai{m}*T{m-1}+Bi{m}*T{m})\(fi{m}-Ai{m}*s{m-1}-Bi{m}*s{m});
for k=2:m
    x{k}=s{k}+T{k}*x{1};
end
    