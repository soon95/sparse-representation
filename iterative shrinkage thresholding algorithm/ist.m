function [ x ] = ist( y,A,lamda,maxErr,maxIter )
%IST 迭代阈值收缩算法
%   输入参数：
%     y:输入信号
%     A:字典
%     lamda:乘子
%     err:收敛误差
%     maxIter:最大迭代次数

    if nargin < 5    
        maxIter = 10000;    
    end    
    if nargin < 4      
        maxErr = 1e-2;      
    end     
    if nargin < 3      
        lamda = 0.1*max(abs(A'*y));     
    end     
    [y_rows,y_columns] = size(y);      
    if y_rows<y_columns      
        y = y';%y should be a column vector      
    end 
    soft = @(x,T) sign(x).*max(abs(x) - T,0);
    n = size(A,2);    
    x = zeros(n,1);%Initialize x=0  
    f = 0.5*(y-A*x)'*(y-A*x)+lamda*sum(abs(x));%added in v1.1
    iter = 0; 
    while 1    
        x_pre = x;
        x = soft(x + A'*(y-A*x),lamda);%update x        
        iter = iter + 1;  
        f_pre = f;%added in v1.1
        f = 0.5*(y-A*x)'*(y-A*x)+lamda*sum(abs(x));%added in v1.1
        if abs(f-f_pre)/f_pre<maxErr%modified in v1.1
            fprintf('abs(f-f_pre)/f_pre<%f\n',maxErr);
            fprintf('IST loop is %d\n',iter);
            break;
        end
        if iter >= maxIter
            fprintf('loop > %d\n',maxIter);
            break;
        end
        if norm(x-x_pre)<maxErr
            fprintf('norm(x-x_pre)<%f\n',maxErr);
            fprintf('IST loop is %d\n',iter);
            break;
        end
    end  
end