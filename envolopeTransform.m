function [ f,p ] = envolopeTransform( sig,fs,isPlot )
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明
% input:
%     sig:输入信号,列向量
%     fs:采样频率
%     isPlot:是否画图，1为是，0为否
%     
% output:
%     f:包络谱频率范围
%     q:各频率分量
%     
% e.g.:
%     [f,p]=envolopeTransform(sig,10240,1);

%%
[sig_rows,sig_columns] = size(sig);  
    if sig_rows<sig_columns  
        sig = sig';            %y should be a column vector  
    end  



[ ~,nfft]=size(sig');
f=(0:nfft/2-1)/nfft*fs;

y=abs(hilbert(sig));
p=abs(fft(y,nfft))/(nfft/2);
p=p(1:nfft/2);


if(isPlot)
    figure();
    plot(f,p);grid on;
    title('包络谱图');
    xlabel('频率 Hz');
    ylabel('幅值 A(m/s^2)'); 
end


end

