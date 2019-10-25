function [ x ] = wist( y,dic,lamda,maxErr,maxIter )
%WIST 迭代阈值收缩算法 软阈值 处理一范数问题
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
        lamda = 0.1*max(abs(dic'*y));     
    end     
    [y_rows,y_columns] = size(y);      
    if y_rows<y_columns      
        y = y';%y should be a column vector      
    end 
    
    % 字典dic为dic_rows*dic_cols矩阵,信号长度为dic_rows,原子个数为dic_cols
    [dic_rows,dic_cols]=size(dic); 
    
    x = zeros(dic_cols,1);%Initialize x=0  
    f = 0.5*(y-dic*x)'*(y-dic*x)+lamda*sum(abs(x));%added in v1.1
    iter = 0; 
    while 1    
        x_pre = x;
        B=x + dic'*(y-dic*x);
        
        % 加权项
        W=weight(y,dic);
        w=W*ones(dic_cols,1);
        
        x = soft_threshold(B,lamda,w);%update x    
        
        %% 查看迭代中的变化情况
%         figure();
%         subplot(3,1,1);
%         plot(x_pre);
%         subplot(3,1,2);
%         plot(x);
%         subplot(3,1,3);
%         plot(B);
        %%
        
        
        iter = iter + 1;  
        f_pre = f;%added in v1.1
        f = 0.5*(y-dic*x)'*(y-dic*x)+lamda*sum(abs(x));%added in v1.1
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

%% 软阈值函数
function [ x ]=soft_threshold(b,lamda,w)
    x=sign(b).*max(abs(b) - lamda*w,0);
end

%% 权值函数
function W=weight(y,dic)
    
    % 字典dic为dic_rows*dic_cols矩阵,信号长度为dic_rows,原子个数为dic_cols
    [dic_rows,dic_cols]=size(dic);
    
    %% 加权方法需要探讨
    % 方案一：采用相关系数，及余弦夹角
    inner_product=1./abs(y'*dic);
    c=inner_product;
    % 方案二：采用峭度指标
    
    
    
    
    % 归一化方法需要探讨
    % c=zscore(inner_product)+1;
    c=mapminmax(c,0,1);
    
    W=diag(c);
    
end
