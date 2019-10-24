function [ x ] = Wbpdn( y,dic,lamda )
%Wbpdn 一种重加权的算法
%   输入参数：
% %     y:输入信号
%     A:字典
%     lamda:乘子

    fprintf('-------开始wbpdn算法-------\n');

    [y_rows,y_columns] = size(y);  
    if y_rows<y_columns  
        y = y';%y should be a column vector  
    end

    % 字典dic为dic_rows*dic_cols矩阵,信号长度为dic_rows,原子个数为dic_cols
    [dic_rows,dic_cols]=size(dic);

    W=weight(y,dic);
    new_w=W*ones(dic_cols,1);
    
    b=dic'*y;
    
    c=lamda*[new_w;new_w]+[-b;b];
    B=[dic'*dic,-dic'*dic;-dic'*dic,dic'*dic];
    lb = zeros(2*dic_cols,1);
    ans=quadprog(B,c,[],[],[],[],lb);
    x=ans(1:dic_cols)-ans(dic_cols+1:2*dic_cols);

end

%% 权值函数
function W=weight(y,dic)
    
    % 字典dic为dic_rows*dic_cols矩阵,信号长度为dic_rows,原子个数为dic_cols
    [dic_rows,dic_cols]=size(dic);

    inner_product=abs(y'*dic);
    
    % 归一化方法需要探讨
%     c=zscore(inner_product)+1;
    c=mapminmax(inner_product,0,1);
    
    W=diag(c);
    
end


