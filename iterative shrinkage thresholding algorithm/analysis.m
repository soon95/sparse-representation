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
A2=1;                 %幅值

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
SNR=-5;


% figure()
% subplot(2,1,1);
% plot(t,sig);
% title('无噪仿真信号');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% ylim([-1,1]);
% 
% subplot(2,1,2);
% plot(t,signal);
% title('无噪仿真信号');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% ylim([-1,1]);

%% 构造字典
f_min=299;                  %(需要根据实际情况调整)
f_max=301;                  %(需要根据实际情况调整)
zeta_min=0.049;              %(需要根据实际情况调整)
zeta_max=0.051;              %(需要根据实际情况调整)
W_step=4;
[Dic,rows,cols]=generate_dic(len,f_min,f_max,zeta_min,zeta_max,W_step,fs);
Dic=dictnormalize(Dic);
Dic=Dic/norm(Dic); 



%% 稀疏恢复算法


%% 比较多种算法重构效果图

% SNR=-6;
% [signal,noise]=noisegen(sig,SNR);
% figure();
% subplot(3,1,1);
% plot(t,sig);
% title('无噪仿真信号');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% ylim([-1,1]);
% subplot(3,1,2);
% plot(t,noise);
% title('噪声');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% ylim([-1,1]);
% subplot(3,1,3);
% plot(t,signal);
% title('仿真信号');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% ylim([-1,1]);
% 
% % CcIST
% sigma=0.02;
% lamda=sigma*sqrt(2*log(cols));
% ts=0.7;
% distance=20;
% theta_CcIST=ClusterShrinkIST(signal,Dic,lamda,distance,ts);
% sig_CcIST=Dic*theta_CcIST;
% 
% % IST
% sigma=0.02;
% lamda=sigma*sqrt(2*log(cols));
% theta_IST=ist(signal,Dic,lamda);
% sig_IST=Dic*theta_IST;
% 
% % CcStOMP
% ts=0.7;
% distance=20;
% maxIter=30;
% [theta_CcStOMP,err_CcStOMP]=ClusterShrinkStOMP(signal,Dic,maxIter,ts,distance);
% sig_CcStOMP=Dic*theta_CcStOMP;
% 
% % 作图
% figure();
% subplot(3,1,1);
% plot(t,sig_CcIST);
% title('CcIST');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% % ylim([-1,1]);
% subplot(3,1,2);
% plot(t,sig_IST);
% title('IST');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% % ylim([-1,1]);
% subplot(3,1,3);
% plot(t,sig_CcStOMP);
% title('CcStOMP');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% % ylim([-1,1]);
%% 比较多种信噪比下相关系数等指标的性能
SNR_range=0:-0.5:-7;
ex_num=50;% 每组实验次数

cc=[];
for SNR=SNR_range
    [signal,noise]=noisegen(sig,SNR);
    
    temp=[];
    for i=1:ex_num
        
        % CcIST
        sigma=0.03;
        lamda=sigma*sqrt(2*log(cols));
        ts=0.7;
        distance=20;
        theta_CcIST=ClusterShrinkIST(signal,Dic,lamda,distance,ts);
        sig_CcIST=Dic*theta_CcIST;

        % IST
        sigma=0.03;
        lamda=sigma*sqrt(2*log(cols));
        theta_IST=ist(signal,Dic,lamda);
        sig_IST=Dic*theta_IST;

        % CcStOMP
        ts=0.7;
        distance=20;
        maxIter=30; 
        [theta_CcStOMP,err_CcStOMP]=ClusterShrinkStOMP(signal,Dic,maxIter,ts,distance);
        sig_CcStOMP=Dic*theta_CcStOMP;
        
        % 算相关度
        r_CcIST=corrcoef(sig,sig_CcIST);
        r_IST=corrcoef(sig,sig_IST);
        r_CcStOMP=corrcoef(sig,sig_CcStOMP);
        
        temp=[temp;r_CcIST(1,2) r_IST(1,2) r_CcStOMP(1,2)];
        
    end
    
    cc=[cc;mean(temp)];
    

end

figure()
plot(SNR_range,cc(:,1),'xr-');
hold on;
plot(SNR_range,cc(:,2),'og-');
plot(SNR_range,cc(:,3),'*b-');
legend('CcIST','IST','CcStOMP');

%%



