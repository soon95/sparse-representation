 clc;clear all;close all;
%% StOMP算法性能 仿真信号分析

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

% figure();
% plot(t,sig);

%% 加入噪声
SNR=-10;
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

%% 构造字典
f_min=299;                  %(需要根据实际情况调整)
f_max=301;                  %(需要根据实际情况调整)
zeta_min=0.049;              %(需要根据实际情况调整)
zeta_max=0.051;              %(需要根据实际情况调整)
W_step=4;
[Dic,rows,cols]=generate_dic(len,f_min,f_max,zeta_min,zeta_max,W_step,fs);
Dic=dictnormalize(Dic);
%% 在进行阈值迭代收缩算法时，需要将字典的矩阵二范数归一化成1；而在其他算法中通常的做法是原子归一化，虽然这样做只是为了好看，但是在阈值迭代收缩算法中由于使用到了majorization-minimization框架，因此对字典尺度有硬性的要求。
% 这个问题困扰了我好久，是重推公式的时候才发现的
Dic=Dic/norm(Dic); 
%% 稀疏恢复 

maxErr=1e-3;
maxIter=100;
window=100;

sigma = 0.025;
lamda = sigma*sqrt(2*log(cols));


[theta_ist,conv_ist]=ist(signal,Dic,lamda,maxErr,maxIter);
sig_recovery_ist=Dic*theta_ist;


[theta_sist,conv_list]=sist(signal,Dic,lamda,maxErr,maxIter,window);
sig_recovery_sist=Dic*theta_sist;


figure()
subplot(4,1,1);
plot(t,sig);
title('signal without noise');
xlabel('Time(s)');
ylabel('Amplitude(m/s^2)');
ylim([-1,1]);

subplot(4,1,2);
plot(t,signal);
title('simulation signal');
xlabel('Time(s)');
ylabel('Amplitude(m/s^2)');
% ylim([-1,1]);

subplot(4,1,3);
plot(t,sig_recovery_sist);
title('LIST');
xlabel('Time(s)');
ylabel('Amplitude(m/s^2)');
ylim([-1,1]);

subplot(4,1,4);
plot(t,sig_recovery_ist);
title('IST');
xlabel('Time(s)');
ylabel('Amplitude(m/s^2)');
ylim([-1,1]);


figure()
subplot(2,1,1);
plot(theta_sist);
title('LIST');
ylabel('Amplitude(m/s^2)');
subplot(2,1,2);
plot(theta_ist);
title('IST');
ylabel('Amplitude(m/s^2)');


figure()
plot(conv_list,'xr-');
hold on;
plot(conv_ist,'og-');
title('convergence');
xlabel('iterations');
% ylabel('');
legend('LIST','IST');
