function [Y,NOISE] = noisegen(X,SNR)
%把白噪声叠加到信号上去：
% noisegen add white Gaussian noise to a signal.
% [Y, NOISE] = NOISEGEN(X,SNR) adds white Gaussian NOISE to X.  The SNR is in dB.
%其中X是纯信号，SNR是要求的信噪比，Y是带噪信号，NOISE是叠加在信号上的噪声。
%rng(1);
NOISE=randn(size(X));
NOISE=NOISE-mean(NOISE);
signal_power = 1/length(X)*sum(X.*X);
noise_variance = signal_power / ( 10^(SNR/10) );
NOISE=sqrt(noise_variance)/std(NOISE)*NOISE;
Y=X+NOISE;
