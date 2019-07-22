function [ result ] = kurtosis( x )
%   此处显示详细说明
% This function simply calculates the summed kurtosis of the input
% input:x;
% ouput:result;

result = mean( (sum((x-ones(size(x,1),1)*mean(x)).^4)/(size(x,1)-2))./(std(x).^4) );

end

