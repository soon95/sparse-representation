function [ t,f,coefs ] = time_frequency_transform_sfft( signal,window,noverlap,nfft,fs,isPlot )
%time_frequency_transform_wavelet 小波变换时频图
% input:
%     signal:输入信号
%     window:窗大小
%     noverlap:重叠大小
%     nfft:离散傅里叶变换点数
%     fs:频率
%     isPlot:是否作图,1表示作图
%     
% output:
%     t:时间
%     f:频率
%     coefs:系数


[coefs, f, t, p] = spectrogram(signal,window,noverlap,nfft,fs);

if(isPlot)
    imagesc(t,f,abs(coefs));
    set(gca,'YDir','normal')
    colorbar;
    xlabel('时间 t/s');
    ylabel('频率 f/Hz');
    title('短时傅里叶时频图');
end

end

