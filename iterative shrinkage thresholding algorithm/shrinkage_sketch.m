clc;clear all;close all;
%%
x=(-2:0.01:2);

lamda=1;

y_soft=sign(x).*max(abs(x)-lamda,0);

y_hard=x;
y_hard(abs(y_hard)<sqrt(lamda))=0;

y_2=x./(1+lamda);

figure();
subplot(1,3,1);
plot(x,y_soft,'linewidth',1.5);
xlim([-2,2]);
ylim([-2,2]);
title('\it{l}_1');

subplot(1,3,2);
plot(x,y_hard,'linewidth',1.5);
xlim([-2,2]);
ylim([-2,2]);
title('\it{l}_0');

subplot(1,3,3);
plot(x,y_2,'linewidth',1.5);
xlim([-2,2]);
ylim([-2,2]);
title('\it{l}_2');