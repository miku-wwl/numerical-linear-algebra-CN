% PE·½·¨³ÌÐò
function [x,iter,res,time]=PEa2(A,B,C,f,a,x,tol)
m=size(B,1); S=cell(m,1); T=cell(m-1,1); 
Ii=eye(size(B{1})); iter=0;
mr=norm(cell2mat(f));
S{1}=B{1}; tic
for i=1:m-1
    T{i}=S{i}\C{i}; Gi=(B{i+1}\A{i+1})*(B{i}\C{i});
    S{i+1}=B{i+1}/(Ii+Gi+a*Gi*Gi);
    N{i+1}=A{i+1}*T{i}+S{i+1}-B{i+1};
end
while (iter<100)
    iter=iter+1;
    z{1}=S{1}\f{1};
    for i=2:m
        z{i}=S{i}\(N{i}*x{i}+f{i}-A{i}*z{i-1});
    end
    x{m}=z{m};
    for i=m-1:-1:1
        x{i}=z{i}-T{i}*x{i+1};
    end
    r{1}=f{1}-(B{1}*x{1}+C{1}*x{2});
    for i=2:m-1
        r{i}=f{i}-(A{i}*x{i-1}+B{i}*x{i}+C{i}*x{i+1});
    end
    r{m}=f{m}-(A{m}*x{m-1}+B{m}*x{m});
    res=norm(cell2mat(r))/mr;
    if res<tol, break; end
end
time=toc;
    


