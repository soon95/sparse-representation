%% 实验验证 外圈故障
clc;clear all;close all
%% 实验数据总时长1秒
fs=12000;               % 振动信号采样频率

total_t=0.2;            % 处理时长

total_N=fs*total_t;     % 总采样点数
point_N=1:total_N;      % 采样点

bias_t=4;               % 偏移时间
bias_N=fs*bias_t;       % 偏移点数

t=point_N/fs;           % 时间

%% 构造字典     169 内圈 频率 2809  阻尼 0.08 , 212 2855
f_min=2852;                  %(需要根据实际情况调整)
f_max=2857;                  %(需要根据实际情况调整)
zeta_min=0.08;              %(需要根据实际情况调整)
zeta_max=0.08;              %(需要根据实际情况调整)
W_step=2;
[Dic,rows,cols]=generate_dic(total_N,f_min,f_max,zeta_min,zeta_max,W_step,fs);
Dic=Dic/norm(Dic); 
%% 读取数据
% load('F:\科研\实验数据\CWRU\212.mat');
% 
% original_signal=X212_DE_time(point_N+bias_N);
% 
% % 原始信号幅值归一化
% original_signal=original_signal/abs(max(original_signal));
% 
% % 加点噪声
% amplitude_noise=0.4;
% noise=amplitude_noise*randn(total_N,1);
% original_signal=original_signal+noise;


%% 读取处理过的信号
load('inner_data1.mat');

%%
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
maxIter=200;
window=300;

lamda=0.06;

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
title('(b)');
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
