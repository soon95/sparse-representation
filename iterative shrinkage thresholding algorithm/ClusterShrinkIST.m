function [ theta,err ] = ClusterShrinkIST( y,dic,lamda,distance,ts,max_iter )  
% ClusterShrinkIST 聚类阈值收缩算法
%   输入参数：
%     y:输入信号
%     A:字典
%     lamda:乘子
%     distance:聚类距离尺度
%     ts:原子挑选相关度阈值
%     maxIter:最大迭代次数

    if nargin<6
        max_iter = 10;%S默认值为10  
    end
    if nargin< 5 
        ts = 0.6;%原子相关度阈值
    end
    if nargin < 4  
        distance=5;
    end  
    
    [dic_rows,dic_cols] = size(dic);%传感矩阵A为M*N矩阵
    
    if nargin < 3  
        lamda=0.05*sqrt(2*log(dic_cols));
    end
    
    [y_rows,y_columns] = size(y);  
    if y_rows<y_columns  
        y = y';            %y should be a column vector  
    end
    
    theta = zeros(dic_cols,1);%用来存储恢复的theta(列向量)  
    Pos_theta = [];%用来迭代过程中存储A被选择的列序号  
    r_n = y;%初始化残差(residual)为y  
    
    err=[r_n];
    
    for ss=1:max_iter%最多迭代S次  
        product = dic'*r_n;%传感矩阵A各列与残差的内积  
        sigma = norm(r_n)/sqrt(dic_rows);%参见参考文献第3页Remarks(3)  
        Js = find(abs(product)>ts*sigma);%选出大于阈值的列  
        
        %% 调试用 对选出的原子进行二次筛选
        [Js_row,~]=size(Js);
        if Js_row==0
            break;
        end
        
        [class,class_num,eachClass]=cluster1D(Js,distance);
        
        Js_new=[];
        for i=1:class_num
            selected=[];
            c=class(i,1:eachClass(i));
            selected(c)=product(c);
            [~,index]=max(abs(selected));
            Js_new(i,1)=index;
        end
        
        %%
        Is = union(Pos_theta,Js_new);%Pos_theta与Js并集  
        if length(Pos_theta) == length(Is)  
            if ss==1  
                theta_ls = 0;%防止第1次就跳出导致theta_ls无定义  
            end  
            break;%如果没有新的列被选中则跳出循环  
        end  
        %At的行数要大于列数，此为最小二乘的基础(列线性无关)  
        if length(Is)<=dic_rows  
            Pos_theta = Is;%更新列序号集合  
            At = dic(:,Pos_theta);%将A的这几列组成矩阵At  
        else%At的列数大于行数，列必为线性相关的,At'*At将不可逆  
            if ss==1  
                theta_ls = 0;%防止第1次就跳出导致theta_ls无定义  
            end  
            break;%跳出for循环  
        end  
        %y=At*theta，以下求theta的最小二乘解(Least Square)  

        
        % CcStOMP中用的是最小二乘解
        % theta_ls = (At'*At)\(At'*y);%最小二乘解  
        % 这里可以将最小二乘解优化为阈值迭代收缩算法
        theta_ls=ist(y,At,lamda);
        

        %At*theta_ls是y在At列空间上的正交投影  
        pre_r_n=r_n;
        r_n = y - At*theta_ls;%更新残差  
        
        err=[err r_n];
        
        if norm(r_n-pre_r_n)/length(r_n)<5e-4%Repeat the steps until r=0  
            break;%跳出for循环  
        end 
    end  
    theta(Pos_theta)=theta_ls;%恢复出的theta  
end