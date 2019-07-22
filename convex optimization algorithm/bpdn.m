function [ x ] = bpdn( y,A,lamda )
%BPDN_quadprog(Basis Pursuit DeNoising with quadprog) Summary of this function goes here
%Version: 1.0 written by jbb0523 @2016-07-22 
%Reference:Nowak R D, Wright S J. Gradient projection for sparse reconstruction:
%Application to compressed sensing and other inverse problems[J]. 
%IEEE Journal of selected topics in signal processing, 2007, 1(4): 586-597.
%(Available at: http://pages.cs.wisc.edu/~swright/papers/FigNW07a.pdf)
    [y_rows,y_columns] = size(y);  
    if y_rows<y_columns  
        y = y';%y should be a column vector  
    end 
    n = size(A,2);
    %according to section II-A of the reference
    b = A'*y;
    c = lamda*ones(2*n,1)+[-b;b];
    B = [A'*A,-A'*A;-A'*A,A'*A];
    lb = zeros(2*n,1);
    z0 = quadprog(B,c,[],[],[],[],lb);
    x = z0(1:n) - z0(n+1:2*n);
end