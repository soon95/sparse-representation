function [ x ] = bpdn( y,A,lamda )
%bpdn 迭代阈值收缩算法
%   输入参数：
%     y:输入信号
%     A:字典
%     lamda:乘子
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
    z0 = quadprog(B,c,[],[],[],[],lb); % 中间四个参数分别是不等于约束和等于约束：A*x<b,A'*x=b'的参数，最后两个参数是参数范围约束，这里要求结果都不小于零
    x = z0(1:n) - z0(n+1:2*n);
end