clc;clear all;close all;
%% 三种范数的图
x=(-2:0.01:2);

y_norm1=abs(x);
y_norm2=x.^2;
y_norm5=abs(x).^0.5;

y_norm0=abs(x).^0.000000000000000000000000000000000000001;

figure();
plot(x,y_norm1,'linewidth',1.5);
hold on;
plot(x,y_norm2,'linewidth',1.5);
plot(x,y_norm5,'linewidth',1.5);
plot(x,y_norm0,'linewidth',1.5);
xlim([-2,2]);
ylim([0,2]);
legend('\it{l}_{1}^{1}','\it{l}_{2}^{2}','\it{l}_{1/2}^{1/2}','\it{l}_{0}^{0}');