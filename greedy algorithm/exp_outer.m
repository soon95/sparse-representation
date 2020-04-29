clc;clear all;close all;
%% 实验研究 参数设置调试程序

% load '深沟球滚动体1加载通道4.mat';
% load '深沟球外圈1x0.5加载通道4.mat';%
% load '深沟球外圈0.2x0.2加载通道4.mat';%从0秒开始，f_min=2245;f_max=2255;zeta_min=0.135;zeta_max=0.138;ts=2.33;  
fs=10240;   %采样频率，手动设置
total_t=0.5;   %采样时间   可以自己设定
total_N=fs*total_t;   %总采样点数
point_N=1:total_N;   %采样点数
bias_t=fs*0;    % 时间偏移   10, 10.5
t=point_N/fs;   %时间

% signal=Data(point_N+bias_t);
%% 滤波



% figure();
% plot(signal);
% fouriorTransform( signal,fs,1);


% bpFilt = designfilt('bandpassfir','FilterOrder',30, ...
%          'CutoffFrequency1',100,'CutoffFrequency2',2800, ...
%          'SampleRate',10240);
% filtered_signal=filter(bpFilt,signal);


% filtered_signal=signal;

% figure();
% plot(filtered_signal);
% fouriorTransform( filtered_signal,fs,1);






%%

length=total_N;

% A_noise=7.2;
% noise=A_noise*randn(length,1);
% filtered_signal=filtered_signal+noise;

load 'outer_data3.mat';  % 3,

%% 构造字典
%外圈信号，原子震荡频率2250，阻尼比0.12
%内圈信号，原子震荡频率1272，阻尼比0.17

f_min=2245;          %内圈，原子震荡频率1272
f_max=2255;
zeta_min=0.119;
zeta_max=0.122;
W_step=1;
[Dic,rows,cols]=generate_dic(length,f_min,f_max,zeta_min,zeta_max,W_step,fs);
Dic=dictnormalize(Dic);
%% 分段进行信号重构

% sig_recover=[];
% ceo_recover=[];

maxIter=100;           %迭代次数
ts=2.7;               %极限系数
distance=100;         %聚类距离尺度，需要着重设置
[theta1]=ClusterShrinkStOMP(filtered_signal,Dic,maxIter,ts,distance);
sig_recover1=Dic*theta1;

[theta2]=StOMP(filtered_signal,Dic,maxIter,ts);
sig_recover2=Dic*theta2;


%%

% figure();
% plot(t,filtered_signal);
% 
% [f0,p0]=envolopeTransform( filtered_signal,fs,0 );
% figure();
% plot(f0,p0)
% xlim([0,1000]);
% ylim([0,0.7]);
% 
% figure();
% plot(t,sig_recover1);
% 
% [f1,p1]=envolopeTransform( sig_recover1,fs,0 );
% figure();
% plot(f1,p1)
% xlim([0,1000]);
% % ylim([0,0.7]);
% 
% 
% figure();
% plot(f0,p0,'b');
% hold on;
% plot(f1,p1,'r')
% xlim([0,1000]);
% % ylim([0,0.7]);

%%
figure();
subplot(3,1,1);

plot(t,filtered_signal)
title('(a)');
xlabel('Time(s)');
ylabel('Amplitude');


[fre2,q2]=fouriorTransform(filtered_signal,fs,0);

subplot(3,1,2);
plot(fre2,q2);
title('(b)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');


[f0,p0]=envolopeTransform( filtered_signal,fs,0 );
subplot(3,1,3);
plot(f0,p0);
title('(c)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
xlim([0,1000]);
ylim([0,1]);




figure();
subplot(3,1,1);
plot(t,sig_recover1);
title('(a)');
xlabel('Time(s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(theta1);
title('(b)');

[f1,p1]=envolopeTransform( sig_recover1,fs,0 );
subplot(3,1,3);
plot(f1,p1);
title('(c)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
xlim([0,1000]);
ylim([0,1]);


figure();
subplot(3,1,1);
plot(t,sig_recover2);
title('(a)');
xlabel('Time(s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(theta2);
title('(b)');

[f2,p2]=envolopeTransform( sig_recover2,fs,0 );
subplot(3,1,3);
plot(f2,p2);
title('(c)');
xlabel('Frequency(Hz)');
ylabel('Amplitude');
xlim([0,1000]);
ylim([0,1]);