clc;clear all;close all;
%% 测试用

fs=2000;                %采样频率
total_time=1;           %采样时间
len=fs*total_time;   %采样点数
point=1:len;         %采样点序列
t=point/fs;             %时间

damp=0.05;              %阻尼系数
f_vibrate=300;          %震荡频率
t0=0.05;                %时延
T=0.1;                  %信号周期


while (t0-T>0)
    t0=t0-T;
end
sig_laplacewavelet=exp(-(damp/sqrt(1-damp^2))*2*pi*f_vibrate*t).*sin(2*pi*f_vibrate*t);
Wss=round(t0*fs);

A1=1;                   %幅值
A2=0.5;                 %幅值

sig=zeros(len,1);
for k=1:len
    if k<=Wss
        sig(k)=0;
    else
        sig(k)=A1*sig_laplacewavelet(k-Wss);
    end
end

count=0;
for i=1:len/(T*fs)
    Wss=round(T*fs*i+t0*fs);
    for k=1:len      
        if k>Wss
            if mod(count,2)==0
                sig(k)=sig(k)+A2*sig_laplacewavelet(k-Wss);
            else
                sig(k)=sig(k)+A1*sig_laplacewavelet(k-Wss);
            end
        end
    end
    count=count+1;
end

%% 加入噪声
SNR=-3;
[signal,noise]=noisegen(sig,SNR);

figure()
subplot(2,1,1);
plot(t,sig);
title('无噪仿真信号');
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');
ylim([-1,1]);

subplot(2,1,2);
plot(t,signal);
title('仿真信号');
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');
ylim([-1,1]);

%%
fouriorTransform(signal,fs,1);

%% 短时傅里叶
figure()
WINDOW = 128;
NOVERLAP = 120;
NFFT = 256;
% MINDB = 20;
% [B, F, T, P] = spectrogram(signal,WINDOW,NOVERLAP,NFFT,fs);
% imagesc(T,F,abs(B));
% set(gca,'YDir','normal')
% colorbar;
% xlabel('时间 t/s');
% ylabel('频率 f/Hz');
% title('短时傅里叶时频图');

time_frequency_transform_sfft(signal,WINDOW,NOVERLAP,NFFT,fs,1);


%% 小波 时频图
% totalscal=256;
% wavename='cmor3-3';
% Fc=centfrq(wavename); % 小波的中心频率  测得Fc = 3
% c=2*Fc*totalscal;    % 测得c = 1536
% scals=c./(1:totalscal);
% f=scal2frq(scals,wavename,1/fs); % 将尺度转换为频率
% coefs = cwt(signal,scals,wavename); % 求连续小波系数
% t=0:1/fs:size(signal)/fs;
% figure()
% imagesc(t,f,abs(coefs));
% set(gca,'YDir','normal')
% colorbar;
% xlabel('时间 t/s');
% ylabel('频率 f/Hz');
% title('小波时频图');

time_frequency_transform_wavelet(signal,256,fs,1);


