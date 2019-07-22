%%
clc;clear all;close all
%————————————————————————————————————————————————
fs=10240;                             %振动信号采样频率(修改)
dt=1/fs;
load '深沟球外圈1x0.5加载通道4.mat';
sig_x=Data(1:1024*2);                          %取信号
N=length(sig_x);
t=1:N;
t1=t/fs;
Ws=200;           %Laplace小波支撑长度（以点数表示）
%%
%进行FFT变换并做频谱图
mag=abs(fft(sig_x));%进行fft变换
figure(1)
plot(t1,sig_x);
title('轴承故障信号');
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');
figure(2)
plot((0:N/2-1)/N*fs,mag(1:N/2));%绘制频谱图
xlabel('频率 F/(Hz)');ylabel('幅值 A/(m/s^2)');
title('信号频谱(FFT)')
%%
display_step=1 %程序执行的现实标志

f_kr=1900:2500;   %小波原子的频率范围
zeta_kr=[0.12:0.005:0.13];     

t_kr_all=N/fs;
t_kr_step=t_kr_all/N;  
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
        sig1=exp(-(zeta/sqrt(1-zeta^2))*2*pi*f*(t/fs));
        for T0=t_kr
            k=k+1;
            Wss=round(T0*fs);
            for m=1:N           %构造Laplace小波原子，延时T0
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
            sig=sig2.*sin(2*pi*f*(t/fs));
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
sig11=exp(-(maxkr_zeta/sqrt(1-maxkr_zeta^2))*2*pi*maxkr_f*(t/fs)).*sin(2*pi*maxkr_f*(t/fs)); %找到的基底
 Wss=round(maxkr_T0*fs);
 for m=1:N
    if m<=Wss
        sig_basis(m)=0;
    else
        sig_basis(m)=sig11(m-Wss); %带延时的基底
    end
 end
figure(3)
plot(t1,sig_basis)
title('basis');
xlabel('time ')
ylabel('amplitude')