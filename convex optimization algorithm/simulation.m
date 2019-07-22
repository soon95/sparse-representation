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

figure();
plot(t,sig);

%% 加入噪声
% A_noise=0.3;
% noise=A_noise*randn(length,1);
% signal=sig+noise;

SNR=-3;
[signal,noise]=noisegen(sig,SNR);


% figure()
% plot(t,sig);
% title('无噪仿真信号');
% xlabel('时间 t/s');
% ylabel('幅值 A(m/s^2)');
% ylim([-1,1]);

% figure();
% plot(t,sig);

%% 构造字典
f_min=299;                  %(需要根据实际情况调整)
f_max=301;                  %(需要根据实际情况调整)
zeta_min=0.049;              %(需要根据实际情况调整)
zeta_max=0.051;              %(需要根据实际情况调整)
W_step=4;
[Dic,rows,cols]=dic(len,f_min,f_max,zeta_min,zeta_max,W_step,fs);
Dic=dictnormalize(Dic);
%% 稀疏恢复

maxIter=30;           %迭代次数
ts=3;               %极限系数
distance=5;         %聚类距离尺度，需要着重设置

%CcStOMP
[theta1,err1]=ClusterShrinkStOMP(signal,Dic,maxIter,ts,distance);
sig_recover1=Dic*theta1;
%StOMP
[theta2,err2]=StOMP(signal,Dic,maxIter,ts);
sig_recover2=Dic*theta2;
%OMP
col=size(Dic,2);
[W,Gamma,err3] = OMP(signal',Dic,maxIter);
theta3=zeros(1,col);
theta3(Gamma)=W;
theta3=theta3';
sig_recover3=Dic*theta3;  


% figure();
% plot(t,sig_recover);
%% 输出图像
figure()
subplot(5,1,1);
plot(t,sig);
title('无噪仿真信号');
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');
ylim([-1,1]);

subplot(5,1,2);
plot(t,signal);
title('含噪仿真信号');
ylim([-2,2]);
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');

subplot(5,1,3);
plot(t,sig_recover1);
title('重构信号');
ylim([-1,1]);
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');

subplot(5,1,4);
plot(t,sig_recover2);
title('重构信号');
ylim([-1,1]);
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');
subplot(5,1,5);
plot(t,sig_recover3);
title('重构信号');
ylim([-1,1]);
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');

figure();
subplot(3,1,1);
plot(theta1);
title('稀疏表示系数');
%ylim([-2,2]);
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');
subplot(3,1,2)
plot(theta2);
title('稀疏表示系数');
%ylim([-2,2]);
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');
subplot(3,1,3)
plot(theta3);
title('稀疏表示系数');
%ylim([-2,2]);
xlabel('时间 t/s');
ylabel('幅值 A(m/s^2)');

%% 相关度
R1=corrcoef(sig,sig_recover1);
display(R1);
R2=corrcoef(sig,sig_recover2);
display(R2);
R3=corrcoef(sig,sig_recover3);
display(R3);
%% 稀疏度
[ceo_num,~]=size(theta1);
sp1=0;
for i=1:ceo_num
    if theta1(i)~=0
        sp1=sp1+1;
    end
end
sp2=0;
for i=1:ceo_num
    if theta2(i)~=0
        sp2=sp2+1;
    end
end
sp3=0;
for i=1:ceo_num
    if theta3(i)~=0
        sp3=sp3+1;
    end
end
display(sp1);
display(sp2);
display(sp3);
%% 残差收敛情况
[max_l,~]=max([size(err1,2),size(err2,2),size(err3,2)]);
for i=1:max_l
   if i<=size(err1,2)
       err(1,i)=norm(err1(:,i));
   else
       err(1,i)=err(1,i-1);
   end
end

for i=1:max_l
   if i<=size(err2,2)
       err(2,i)=norm(err2(:,i));
   else
       err(2,i)=err(2,i-1);
   end
end

for i=1:max_l
   if i<=size(err3,2)
       err(3,i)=norm(err3(:,i));
   else
       err(3,i)=err(3,i-1);
   end
end

err_noise=ones(1,max_l).*norm(noise);

err=[err;err_noise];

figure();
subplot(2,1,1);
plot([0:max_l-1],err(1,:),'xr-');
hold on;
plot([0:max_l-1],err(2,:),'og-');
plot([0:max_l-1],err(3,:),'*b-');
plot([0:max_l-1],err(4,:),'k-');
title('(a)');
xlabel('Iteration');
ylabel('E');
legend('CcStOMP','StOMP','OMP','noise');
%% 冲击成分提取情况，峭度值
[max_l,~]=max([size(err1,2),size(err2,2),size(err3,2)]);

for i=1:max_l
   if i<=size(err1,2)
       kur(1,i)=kurtosis(err1(:,i));
   else
       kur(1,i)=kur(1,i-1);
   end
end

for i=1:max_l
   if i<=size(err2,2)
       kur(2,i)=kurtosis(err2(:,i));
   else
       kur(2,i)=kur(2,i-1);
   end
end

for i=1:max_l
   if i<=size(err3,2)
       kur(3,i)=kurtosis(err3(:,i));
   else
       kur(3,i)=kur(3,i-1);
   end
end

kur_noise=ones(1,max_l).*kurtosis(noise);

kur=[kur;kur_noise];

subplot(2,1,2);
plot([0:max_l-1],kur(1,:),'xr-');
hold on;
plot([0:max_l-1],kur(2,:),'og-');
plot([0:max_l-1],kur(3,:),'*b-');
plot([0:max_l-1],kur(4,:),'k-');
title('(b)');
xlabel('Iteration');
ylabel('K');
legend('CcStOMP','StOMP','OMP','noise');