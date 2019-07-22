function [ f,q ] = fouriorTransform( sig,fs,isPlot)
%fouriorTransform 此处显示有关此函数的摘要
%   说明：对信号进行快速傅里叶变换，并作图
% input:
%     sig:输入信号
%     fs:采样频率
%     isPlot:是否画图，1为是，0为否
%     
% output:
%     f:频谱图频率范围
%     q:各频率分量
%     
% e.g.:
%     [f,q]=fouriorTransform(sig,10240,1);


%%
[ ~,nfft]=size(sig');
f=(0:nfft/2-1)/nfft*fs;
q=abs(fft(sig,nfft))/(nfft/2);
q=q(1:nfft/2);

if(isPlot)
    figure();
    plot(f,q(1:nfft/2));grid on;
    title('原始故障频域图');
    xlabel('频率 Hz');
    ylabel('幅值 A(m/s^2)');
end

end

