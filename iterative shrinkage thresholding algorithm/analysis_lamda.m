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
[Dic,rows,cols]=generate_dic(length,f_min,f_max,zeta_min,zeta_max,W_step,fs);
Dic=dictnormalize(Dic);
Dic=Dic/norm(Dic); 
%% 
% SNR=-9;
ex_num=50;          %实验次数  50

maxErr=1e-3;        %误差限
maxIter=100;        %最大迭代次数
window=100;         %窗口大小

% sigma = 0.025;
% lamda = sigma*sqrt(2*log(cols));


snr_inter=-3:-1:-14;

lamda_range=0.05:0.05:0.25;

change_list=[];
change_ist=[];

for lamda=lamda_range

    fin=[];  

    for SNR=snr_inter

        CC1=[];
        CC2=[];

        for i=1:ex_num
            [signal,~]=noisegen(sig,SNR);

            %list
            [theta_list,~]=sist(signal,Dic,lamda,maxErr,maxIter,window);
            sig_recovery_list=Dic*theta_list;
            %ist
            [theta_ist,~]=ist(signal,Dic,lamda,maxErr,maxIter);
            sig_recovery_ist=Dic*theta_ist;



            %% 性能统计
            % 相似性
            r1=corrcoef(sig,sig_recovery_list);
            CC1=[CC1 r1(1,2)];
            r2=corrcoef(sig,sig_recovery_ist);
            CC2=[CC2 r2(1,2)];    


        end

        res=[];
        res(1,:)=CC1;
        res(2,:)=CC2;


        ans=mean(res,2);
        % display(ans);



        fin=[fin ans];

    end

    change_list=[change_list fin(1,:)'];
    change_ist=[change_ist fin(2,:)'];


end
%%
upper_list=max(change_list,[],2);
lower_list=min(change_list,[],2);
change_list=[change_list upper_list lower_list];

upper_ist=max(change_ist,[],2);
lower_ist=min(change_ist,[],2);
change_ist=[change_ist upper_ist lower_ist];

figure()
subplot(1,2,1);
plot(snr_inter,change_list);
title('(a)');
xlabel('SNR/dB');
set(gca, 'XDir','reverse'); 
ylabel('ACC');
legend('lamda=0.05','lamda=0.1','lamda=0.15','lamda=0.2','lamda=0.25','upper bound','lower bound');


subplot(1,2,2);
plot(snr_inter,change_ist);
title('(b)');
xlabel('SNR/dB');
set(gca, 'XDir','reverse'); 
ylabel('ACC');
legend('lamda=0.05','lamda=0.1','lamda=0.15','lamda=0.2','lamda=0.25','upper bound','lower bound');

