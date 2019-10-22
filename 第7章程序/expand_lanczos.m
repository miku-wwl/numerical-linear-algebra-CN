function [V,alpha,beta]=expand_lanczos(A,V,alpha,beta,m,k)
%�����Գƾ���A�ĳ���Ϊm��Lanczos�ֽ�,���㷨������չΪһ������Ϊk��Lanczos�ֽ�
for i=m+1:k
    u=A*V(:,i)-beta(i-1)*V(:,i-1);
    alpha(i)=u'*V(:,i);
    u=u-alpha(i)*V(:,i);
    u=u-(u'*V(:,i-1))*V(:,i-1); %�ֲ���������
    u=u-(u'*V(:,i))*V(:,i);
    u=u-V(:,1:i)*(V(:,1:i)'*u); %��ȫ��������
    beta(i)=norm(u);
    if (beta(i)==0)
        return;
    else
        V(:,i+1)=u/beta(i);
        V=[V, V(:,i+1)];
    end
end


    
    


