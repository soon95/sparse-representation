function [A,rows,cols]=dic(N,f_min,f_max,zeta_min,zeta_max,W_step,fs)
% fs=10240;         %采样频率
t=1:N;          %时间长度为10s
rows=length(t); %字典的行数
A1=1;           %有效信号幅值
A=[];
%Wss=round(T*fs);
l1=0;
for Wss=0:W_step:N-1
    l1=l1+1;
    l2=0;
   for  zeta0=zeta_min:.01:zeta_max     %(需要根据实际情况调整)
       l2=l2+1;
       l3=0;
       for f0=f_min:f_max
           l3=l3+1;
         sig1=exp(-(zeta0/sqrt(1-zeta0^2))*2*pi*f0*(t/fs)).*sin(2*pi*f0*(t/fs));
          for k=1:N
               if k<=Wss
                   sig2(k)=0;
               else
                   sig2(k)=sig1(k-Wss);
               end
          end
            atom=sig2;
            A=[A atom'];
       end
   end
end
cols=l1*l2*l3;         %字典的列数
end
