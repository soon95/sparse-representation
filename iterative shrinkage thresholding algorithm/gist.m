function [ theta , err ] = gist( y,dic,lamda,max_reco_iter,max_err,ts,max_sel_iter )
%GIST 贪婪迭代收缩算法
%   输入参数：
%     y:输入信号
%     A:字典
%     lamda:乘子
%     err:收敛误差
%     maxIter:最大迭代次数
%   输出参数：
%     x:稀疏系数向量
%     err:残差随迭代变化情况

    fprintf('-------开始gist算法-------\n');

    if nargin<7
        max_sel_iter=10;% 原子选择中的迭代次数
    end
    
    if nargin<6
        ts=1;%ts范围[2,3],默认值为2.5 
    end

    if nargin<5
        max_err=1e-2;% 收敛性判断
    end
    
    if nargin<4
       max_reco_iter=10000;% 阈值迭代收缩中的迭代次数
    end
    
    if nargin<3
        lamda=0.01;% 拉格朗日乘子
    end
    
    [y_rows,y_cols]=size(y);
    % 确保 y 是一个列向量
    if y_rows<y_cols
        y=y';
    end
    
    % 字典dic为dic_rows*dic_cols矩阵,信号长度为dic_rows,原子个数为dic_cols
    [dic_rows,dic_cols]=size(dic);
    
    % 初始化theta向量
    theta=zeros(dic_cols,1);
    
    % 储存被选中原子列序号
    index_theta=[];
    
    % 初始化残差
    res=y;
    err=[res];
    
    % 最多迭代max_sel_iter次
    for sel_iter=1:max_sel_iter
        product=dic'*res;% 计算字典中原子与残差的内积
        
        % 挑选出大于阈值的列
        sigma=norm(res)/sqrt(dic_rows);
        selected=find(abs(product)>ts*sigma);
        % 通常的迭代停止条件在这里，当没有挑选出新的原子时停止迭代
        if isempty(selected)
            fprintf('sel_iter:%f-更新支撑集为空，跳出循环，完成重构\n',sel_iter);
            break;
        end
        index_theta=union(index_theta,selected);
        
        % 生成贪婪字典，用此字典进行迭代阈值收缩
        greedy_dic=dic(:,index_theta);
        [~,greedy_dic_cols]=size(greedy_dic);
        sel_theta=zeros(greedy_dic_cols,1);
        
        % 目标函数收敛性
        obj_f=0.5*(y-greedy_dic*sel_theta)'*(y-greedy_dic*sel_theta)+lamda*sum(abs(sel_theta));
        
        for reco_iter=1:max_reco_iter
            sel_theta_pre=sel_theta;
            
            % 系数更新阶段
            sel_theta=shrinkage_function(sel_theta+greedy_dic'*(y-greedy_dic*sel_theta),lamda);
            
            obj_f_pre=obj_f;
            obj_f=0.5*(y-greedy_dic*sel_theta)'*(y-greedy_dic*sel_theta)+lamda*sum(abs(sel_theta));
            
            if abs((obj_f-obj_f_pre)/obj_f_pre)<max_err
                fprintf('sel_iter:%f-abs((obj_f-obj_f_pre)/obj_f_pre)<%f\n',sel_iter,max_err);
                break; 
            end
            
            if norm(sel_theta-sel_theta_pre)<max_err
                fprintf('sel_iter:%f-norm(x-x_pre)<%f\n',sel_iter,max_err);
                break;
            end
            
            % 这里可能需要补充其他跳出规则
            
        end
        
        res_pre=res;
        res=y-greedy_dic*sel_theta;
        
        err=[err res];
        
        if abs((norm(res)-norm(res_pre))/norm(res_pre))<max_err
            fprintf('sel_iter:%f-残差收敛，跳出循环，完成重构\n',sel_iter);
            break;
        end    
    end
    
    % 恢复出的theta
    theta(index_theta)=sel_theta;

end

%% 收缩函数
function [x] =shrinkage_function(b,lamda)
     x=sign(b).*max(abs(b) - lamda,0);
end

