clc;clear all;close all;
%% 实验研究 参数设置调试程序

% load '深沟球滚动体1加载通道4.mat';
% load '深沟球外圈0.2x0.2加载通道4.mat';%从0秒开始，f_min=2245;f_max=2255;zeta_min=0.135;zeta_max=0.138;ts=2.33;  
fs=10240;   %采样频率，手动设置
total_t=0.5;   %采样时间   可以自己设定
total_N=fs*total_t;   %总采样点数
point_N=1:total_N;   %采样点数
bias_t=fs*2.5;    % 时间偏移   10, 10.5
t=point_N/fs;   %时间

%% 构造字典
%外圈信号，原子震荡频率2250，阻尼比0.12
%内圈信号，原子震荡频率1272，阻尼比0.17
% 
% f_min=2245;          %内圈，原子震荡频率1272
% f_max=2255;
% zeta_min=0.119;
% zeta_max=0.122;
% W_step=2;
% [Dic,rows,cols]=generate_dic(total_N,f_min,f_max,zeta_min,zeta_max,W_step,fs);
% Dic=dictnormalize(Dic);
% Dic=Dic/norm(Dic);

% 构造字典太耗时了，把它持久化一下吧
load('Dic_outer2.mat');


%% 读取数据

% load '深沟球外圈1x0.5加载通道4.mat';
% original_signal=Data(point_N+bias_t);
% 
% % 原始信号幅值归一化
% original_signal=original_signal/abs(max(original_signal));
% 
% 
% % 加点随机噪声看看
% amplitude_noise=0.7;
% noise=amplitude_noise*randn(total_N,1);
% original_signal=original_signal+noise;

%% 读取处理过的信号
load 'outer2_data3.mat';
% 
% % 加点随机噪声看看
% amplitude_noise=0.2;
% noise=amplitude_noise*randn(total_N,1);
% original_signal=original_signal+noise;

%% 原始信号

figure();
subplot(3,1,1);

plot(t,original_signal)
title('(a)');
xlabel('Time(s)');
ylabel('Amplitude');


[f1,q_orig]=fouriorTransform(original_signal,fs,0);

subplot(3,1,2);
plot(f1,q_orig);
title('(b)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');


[f2,p_orig]=envolopeTransform( original_signal,fs,0 );
subplot(3,1,3);
plot(f2,p_orig);
title('(c)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
xlim([0,1000]);
ylim([0,0.15]);

%% 重构参数设置

maxErr=1e-4;
maxIter=100;
window=700;     % 这个参数至关重要

lamda=0.12;

%% IST信号重构

theta_ist=ist(original_signal,Dic,lamda,maxErr,maxIter);
sig_recovery_ist=Dic*theta_ist;

figure();
subplot(3,1,1);

plot(t,sig_recovery_ist)
title('(a)');
xlabel('Time(s)');
ylabel('Amplitude');


subplot(3,1,2);
plot(theta_ist);
title('(b)');
xlabel('Index');
ylabel('Amplitude');


[f2,p_ist]=envolopeTransform( sig_recovery_ist,fs,0 );
subplot(3,1,3);
plot(f2,p_ist);
title('(c)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
xlim([0,1000]);
ylim([0,0.1]);


%% LIST信号重构

theta_sist=sist(original_signal,Dic,lamda,maxErr,maxIter,window);
sig_recovery_sist=Dic*theta_sist;

figure();
subplot(3,1,1);

plot(t,sig_recovery_sist)
title('(a)');
xlabel('Time(s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(theta_sist);
title('(b)')
xlabel('Index');
ylabel('Amplitude');


[f2,p_list]=envolopeTransform( sig_recovery_sist,fs,0 );
subplot(3,1,3);
plot(f2,p_list);
title('(c)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
xlim([0,1000]);
ylim([0,0.1]);

