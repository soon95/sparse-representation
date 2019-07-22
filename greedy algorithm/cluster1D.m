function [ class,class_num,eachClass ] = cluster1D( ceo,distance )
%cluster1D 函数摘要
%   
% input:
%     ceo;稀疏系数索引
%     distance:距离，建议设置为5
% output:
%     class:分类情况
%     class_num:总功类数
%     eachClass:每类各有几个


[length,~]=size(ceo);

class=[];
class_num=1;

count=1;
class(class_num,count)=ceo(1,1);
for i=2:length
    temp=ceo(i-1,1);
    cur=ceo(i,1);
    
    if abs(temp-cur)<=distance             %如果相邻，则加入原集合
        count=count+1;
        class(class_num,count)=cur;
        
    else                            %如果不相邻，则创建新集合
        eachClass(class_num)=count;
        
        class_num=class_num+1;
        count=1;
        
        class(class_num,count)=cur;
    end
    
        
end

eachClass(class_num)=count;

end

