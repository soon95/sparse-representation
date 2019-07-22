clc;clear all;close all;
%% StOMP算法性能 仿真信号分析

fs=2000;                %采样频率
total_time=1;           %采样时间
length=fs*total_time;   %采样点数
point=1:length;         %采样点序列
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

sig=zeros(length,1);
for k=1:length
    if k<=Wss
        sig(k)=0;
    else
        sig(k)=A1*sig_laplacewavelet(k-Wss);
    end
end

count=0;
for i=1:length/(T*fs)
    Wss=round(T*fs*i+t0*fs);
    for k=1:length      
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

%% 构造字典
f_min=299;                  %(需要根据实际情况调整)
f_max=301;                  %(需要根据实际情况调整)
zeta_min=0.049;              %(需要根据实际情况调整)
zeta_max=0.051;              %(需要根据实际情况调整)
W_step=4;
[Dic,rows,cols]=dic(length,f_min,f_max,zeta_min,zeta_max,W_step,fs);
Dic=dictnormalize(Dic);

%% 
% A_noise=0.3;
% noise=A_noise*randn(length,1);
% signal=sig+noise;




% SNR=-9;
ex_num=50;          %实验次数  50

maxIter=50;           %迭代次数
ts=3;               %极限系数
distance=5;         %聚类距离尺度，需要着重设置


fin=[];
snr_inter=-3:-0.5:-14;
for SNR=snr_inter

CC1=[];
CC2=[];
CC3=[];
CC4=[];

SP1=[];
SP2=[];
SP3=[];
SP4=[];  
    
for i=1:ex_num
    [signal,~]=noisegen(sig,SNR);
    
    %CcStOMP
    [theta1]=ClusterShrinkStOMP(signal,Dic,maxIter,ts,distance);
    sig_recover1=Dic*theta1;
    %StOMP
    [theta2]=StOMP(signal,Dic,maxIter,ts);
    sig_recover2=Dic*theta2;
    %OMP
    col=size(Dic,2);
    [W,Gamma,error] = OMP(signal',Dic,maxIter);
    theta3=zeros(1,col);
    theta3(Gamma)=W;
    theta3=theta3';
    sig_recover3=Dic*theta3;    
    %CoSaMP
    [theta4]=CoSaMP(signal,Dic,10 );
    sig_recover4=Dic*theta4;
    
    
    %% 性能统计
    r1=corrcoef(sig,sig_recover1);
    CC1=[CC1 r1(1,2)];
    r2=corrcoef(sig,sig_recover2);
    CC2=[CC2 r2(1,2)];    
    r3=corrcoef(sig,sig_recover3);
    CC3=[CC3 r3(1,2)]; 
    r4=corrcoef(sig,sig_recover4);
    CC4=[CC4 r4(1,2)]; 
    
    [ceo_num,~]=size(theta1);
    sp1=0;
    for i=1:ceo_num
        if theta1(i)~=0
            sp1=sp1+1;
        end
    end
    SP1=[SP1 sp1/cols];
    sp2=0;
    for i=1:ceo_num
        if theta2(i)~=0
            sp2=sp2+1;
        end
    end
    SP2=[SP2 sp2/cols];
    sp3=0;
    for i=1:ceo_num
        if theta3(i)~=0
            sp3=sp3+1;
        end
    end  
    SP3=[SP3 sp3/cols];
    
    sp4=0;
    for i=1:ceo_num
        if theta4(i)~=0
            sp4=sp4+1;
        end
    end  
    SP4=[SP4 sp4/cols];
    
end

res=[];
res(1,:)=CC1;
res(2,:)=CC2;
res(3,:)=CC3;
res(4,:)=CC4;
res(5,:)=SP1;
res(6,:)=SP2;
res(7,:)=SP3;
res(8,:)=SP4;

ans=mean(res,2);
% display(ans);


fin=[fin ans];

end

%% 绘图

figure();
subplot(2,1,1);
plot(snr_inter,fin(1,:),'xr-');
hold on;
plot(snr_inter,fin(2,:),'og-');
plot(snr_inter,fin(3,:),'*b-');
plot(snr_inter,fin(4,:),'dc-');
title('(a)');
xlabel('SNR/dB');
set(gca, 'XDir','reverse'); 
ylabel('ACC');
legend('CcStOMP','StOMP','OMP','CoSaMP');


subplot(2,1,2);
plot(snr_inter,fin(5,:),'xr-');
hold on;
plot(snr_inter,fin(6,:),'og-');
plot(snr_inter,fin(7,:),'*b-');
plot(snr_inter,fin(8,:),'dc-');
title('(b)');
xlabel('SNR/dB');
set(gca, 'XDir','reverse'); 
ylabel('ASF');
legend('CcStOMP','StOMP','OMP','CoSaMP');