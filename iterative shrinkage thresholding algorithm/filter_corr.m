%%
clc;clear all;close all
%% 实验数据总时长1秒
fs=12000;               % 振动信号采样频率

total_t=0.1;            % 处理时长

total_N=fs*total_t;     % 总采样点数
point_N=1:total_N;      % 采样点

bias_t=0;               % 偏移时间
bias_N=fs*bias_t;       % 偏移点数

t=point_N/fs;           % 时间
Ws=200;           %Laplace小波支撑长度（以点数表示）


load('F:\科研\实验数据\CWRU\169.mat');

sig_x=X169_DE_time(point_N+bias_N);

%%
%进行FFT变换并做频谱图
mag=abs(fft(sig_x));%进行fft变换
figure(1)
plot(t,sig_x);
title('轴承故障信号'); 
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');
figure(2)
plot((0:total_N/2-1)/total_N*fs,mag(1:total_N/2));%绘制频谱图
xlabel('频率 F/(Hz)');ylabel('幅值 A/(m/s^2)');
title('信号频谱(FFT)')
%%
display_step=1 %程序执行的现实标志

f_kr=2800:2820;   %小波原子的频率范围
zeta_kr=[0.05:0.01:0.09];     

t_kr_all=total_N/fs;
t_kr_step=t_kr_all/total_N;  
t_kr=0:t_kr_step:t_kr_all;
%t_kr=0:t_kr_step:0.015; %%小波原子的时间范围

i=0;j=0;k=0;
maxkr=0;
maxkr_f=min(f_kr);
maxkr_zeta=min(zeta_kr);
maxkr_T0=min(t_kr);
display_step=2.1  %程序执行的现实标志
for f=f_kr
    i=i+1;
    j=0;
    for zeta=zeta_kr
        j=j+1;
        k=0;
        sig1=exp(-(zeta/sqrt(1-zeta^2))*2*pi*f*(point_N/fs));
        for T0=t_kr
            k=k+1;
            Wss=round(T0*fs);
            for m=1:total_N           %构造Laplace小波原子，延时T0
                if m>Wss
                    if m<Wss+Ws
                        sig2(m)=sig1(m-Wss);
                    else 
                        sig2(m)=0;
                    end

                else
                    sig2(m)=0;
                end
            end
            sig=sig2.*sin(2*pi*f*(point_N/fs));
            kr(i,j,k)=abs(dot(sig_x,sig))/(sqrt(sum(sig_x.*sig_x))*sqrt(sum(sig.*sig)));
            if kr(i,j,k)>maxkr
                maxkr=kr(i,j,k);
                maxkr_f=f;          
                maxkr_zeta=zeta;    
                maxkr_T0=T0;                                  
            end     
        end
    end
    display_f_kr=f %程序执行的显示标志
end
%%
display_f=maxkr_f
display_zeta= maxkr_zeta
display_To=maxkr_T0
%% matched basis
sig11=exp(-(maxkr_zeta/sqrt(1-maxkr_zeta^2))*2*pi*maxkr_f*(point_N/fs)).*sin(2*pi*maxkr_f*(point_N/fs)); %找到的基底
 Wss=round(maxkr_T0*fs);
 for m=1:total_N
    if m<=Wss
        sig_basis(m)=0;
    else
        sig_basis(m)=sig11(m-Wss); %带延时的基底
    end
 end
figure(3)
plot(t,sig_basis)
title('basis');
xlabel('time ')
ylabel('amplitude')