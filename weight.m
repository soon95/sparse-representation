function [ W ] = weight( y,dic )
%weight 全指函数
%   




    % 字典dic为dic_rows*dic_cols矩阵,信号长度为dic_rows,原子个数为dic_cols
    [dic_rows,dic_cols]=size(dic);
    
    %% 加权方法需要探讨
    % 方案一：采用相关系数加权
%     c=w_inner_product(y,dic);
    % 方案二：采用峭度指标
    c=w_kurt(y,dic);
    
    
    % 归一化方法需要探讨
    % c=zscore(inner_product)+1;
    % 归一化到01区间
%     c=mapminmax(c,0,1);
    % 均值调为1
    c=c-mean(c)+1;
    
    
    
    W=diag(c);

end

%% 采用相关系数进行加权
function c=w_inner_product(y,dic)

    c=abs(y'*dic);

end

%% 采用峭度进行加权
function c=w_kurt(y,dic)
    
    [dic_rows,dic_cols]=size(dic);
    for i=1:dic_cols
        atom=dic(:,i);
        
        c(i)=kurtsis(y(find(atom~=0)));
        
    end
    
%     c=(1./c).^2;

    

end


