function [ t,f,coefs ] = time_frequency_transform_wavelet( signal,totalScal,fs,isPlot )
%time_frequency_transform_wavelet 小波变换时频图
% input:
%     signal:输入信号
%     totalScal:尺度b
%     fs:采样频率  
%     isPlot:是否作图,1表示作图
%     
% output:
%     t:时间
%     f:频率
%     coefs:系数

wavename='cmor3-3';
Fc=centfrq(wavename); % 小波的中心频率
c=2*Fc*totalScal;
scals=c./(1:totalScal);
f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
coefs = cwt(signal,scals,wavename); % 求连续小波系数
t=0:1/fs:size(signal)/fs;

if(isPlot)
    figure()
    imagesc(t,f,abs(coefs));
    set(gca,'YDir','normal')
    colorbar;
    xlabel('时间 t/s');
    ylabel('频率 f/Hz');
    title('小波时频图');
end

end

