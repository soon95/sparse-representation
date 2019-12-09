function [ x ] = sist( y,A,lamda,maxErr,maxIter,window )
%IST 迭代阈值收缩算法 的一个改进版本
%   输入参数：
%     y:输入信号
%     A:字典
%     lamda:乘子
%     err:收敛误差
%     maxIter:最大迭代次数

    if nargin < 6    
        window = 20;    
    end 
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
    %soft = @(x,T) sign(x).*max(abs(x) - T,0);
    
    n = size(A,2);    
    x = zeros(n,1);%Initialize x=0  
    f = 0.5*(y-A*x)'*(y-A*x)+lamda*sum(abs(x));%added in v1.1
    iter = 0; 
    while 1    
        x_pre = x;
        B=x + A'*(y-A*x);
        x = soft_threshold_v2(B,lamda,window);%update x    
        
        %% 查看迭代中变化情况
        if mod(iter,10)==0
            figure();
            subplot(3,1,1);
            plot(x_pre);
            subplot(3,1,2);
            plot(x);
            subplot(3,1,3);
            plot(B);
        end
        
        %%
        iter = iter + 1;  
        f_pre = f;%added in v1.1
        f = 0.5*(y-A*x)'*(y-A*x)+lamda*sum(abs(x));%added in v1.1
        if abs(f-f_pre)/f_pre<maxErr%modified in v1.1
            fprintf('abs(f-f_pre)/f_pre<%f\n',maxErr);
            fprintf('SIST loop is %d\n',iter);
            break;
        end
        if iter >= maxIter
            fprintf('loop > %d\n',maxIter);
            break;
        end
        if norm(x-x_pre)<maxErr
            fprintf('norm(x-x_pre)<%f\n',maxErr);
            fprintf('SIST loop is %d\n',iter);
            break;
        end
    end  
end

%% 局部软阈值函数
function [ x ]=soft_threshold(b,lamda,window)

    length=size(b,1);
    % 初始化结果，全零
    x=zeros(length,1);
    % 找出大于lamda的项索引
    candi=find(abs(b)>lamda);
    [class,class_num,eachClass]=cluster1D(candi,window);
    
    for i=1:class_num
        atoms_index=class(i,1:eachClass(i));
        atoms=b(atoms_index);
        [~,selected_index]=max(abs(atoms));
        % 在字典中全局的索引位置
        index=atoms_index(selected_index);
        selected=b(index);
        
        val=sign(selected)*(abs(selected)-lamda);
        
        x(index)=val;
        
    end

%     x=sign(b).*max(abs(b) - lamda,0);
end


%% 局部软阈值函数 的另一种实现
function [b]=soft_threshold_v2(b,lamda,window)
    
    length=size(b,1);
    half=round(window/2);
    
    % 首先收缩小于lamda的值
    b(b<lamda)=0;
    
    % 
    residual=abs(b);
    
    for i=1:length
        
        [value,index]=max(residual);
        
        if value==0
            break;
        end
        
        b(index-half:index-1)=0;
        b(index)=sign(b(index))*(value-lamda);
        b(index+1:index+half)=0;
        
        residual(index-half:index+half)=0;
        
    end
end


