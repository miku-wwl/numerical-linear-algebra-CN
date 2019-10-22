%matirx A
function [A,b]=mat(n)
q=ceil(0.9*n);
W=zeros(q); N=zeros(n-q); 
v=zeros(n-q,1); 
for k=1:q
    for j=1:q
        if j==k
            W(k,j)=k+1;
        else if abs(k-j)==1
                W(k,j)=1;
            else
                W(k,j)=0;
            end
        end
    end
end
for k=1:(n-q)
    v(k)=1.0/k;
    for j=1:(n-q)
        if j==k
            N(k,j)=k+1;
        else if abs(k-j)==1
                N(k,j)=1;
            else
                N(k,j)=0;
            end
        end
    end
end
F=zeros(q,n-q); 
for k=1:q
    for j=1:n-q
        if k==j+2*q-n
            F(k,j)=j;
        end
    end
end
O=diag(v);
A=[W, F*O; -F', N];
e=ones(n,1);
b=A*e;

