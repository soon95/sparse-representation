function [ x, obj_f ] = ist( y,A,lamda,maxErr,maxIter )
%IST 迭代阈值收缩算法 软阈值 处理一范数问题
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
        maxErr = 1e-3;      
    end     
    if nargin < 3      
        lamda = 0.1*max(abs(A'*y));     
    end     
    [y_rows,y_columns] = size(y);      
    if y_rows<y_columns      
        y = y';%y should be a column vector      
    end 
    
    obj_f=[];
    n = size(A,2);    
    x = zeros(n,1);%Initialize x=0  
    f = 0.5*(y-A*x)'*(y-A*x)+lamda*sum(abs(x));%added in v1.1
    obj_f=[obj_f f];
    iter = 0; 
    while 1    
        x_pre = x;
        B=x + A'*(y-A*x);
        x = soft_threshold(B,lamda);%update x    
        
        %% 查看迭代中变化情况
%        if mod(iter,10)==0
%             figure();
%             subplot(3,1,1);
%             plot(x_pre);
%             subplot(3,1,2);
%             plot(x);
%             subplot(3,1,3);
%             plot(B);
%         end
        %%
        
        
        iter = iter + 1;  
        f_pre = f;%added in v1.1
        f = 0.5*(y-A*x)'*(y-A*x)+lamda*sum(abs(x));%added in v1.1
        % 更新目标函数值
        obj_f=[obj_f f];
        
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
function [ x ]=soft_threshold(b,lamda)
    x=sign(b).*max(abs(b) - lamda,0);
end
