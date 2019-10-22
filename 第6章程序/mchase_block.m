%��׷�Ϸ�����--mchase_block.m
function [fi]=mchase_block(Ai,Bi,Ci,fi,m)
%�ÿ�׷�Ϸ�������ԽǷ����� Ax=f
%����: AiΪA�Ĵ��¶Խǿ�,BiΪA�����Խǿ�,
%CiΪA�Ĵ��϶Խǿ�, fiΪ�Ҷ�����,mΪ�ֿ���. 
%���: ������fi,����LU�ֽ��L{k},U{k}�����Bi{k},
%Ci{k}��λ��,y{k}��x{k}�Ⱥ�����fi{k}��λ��
fi{1}=Bi{1}\fi{1};
for k=2:m  
    Ci{k-1}=Bi{k-1}\Ci{k-1};
    Bi{k}=Bi{k}-Ai{k}*Ci{k-1};
    fi{k}=Bi{k}\(fi{k}-Ai{k}*fi{k-1});
end  
for k=m-1:-1:1
    fi{k}=fi{k}-Ci{k}*fi{k+1};
end
    