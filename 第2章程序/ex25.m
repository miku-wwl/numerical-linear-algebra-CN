%Àý2.5-ex25
A=[9,1,4,5;6,2,1,3;4,1,9,9;5,2,9,4;4,2,5,1];
[Q1,R1]=G_Schmidt1(A);
err1=norm(A-Q1*R1)
[Q2,R2]=G_Schmidt2(A);
err2=norm(A-Q2*R2)