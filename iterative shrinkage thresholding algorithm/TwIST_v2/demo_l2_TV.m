function demo_l2_TV
% TV based image restoration using TwIST
% cameraman, blur uniform 9*9, BSNR = 40 dB
%
% This demo illustrates the computation of  
%
%     xe = arg min 0.5 ||Ax-y||^2 + tau TV(x)
%             x
% 
% with the TwIST algorithm. 
%function demo_l2_TV
% TV based image restoration using TwIST
% cameraman, blur uniform 9*9, BSNR = 40 dB
%
% This demo illustrates the computation of  
%
%     xe = arg min 0.5 ||Ax-y||^2 + tau TV(x)
%             x
% 
% with the TwIST algorithm. 
%
% The proximal operator, i.e., the solution of
%
%     xe = arg min 0.5 ||x-y||^2 + tau TV(x)
%             x
% is  given the Chambolle algorithm
% 
% A. Chambolle, “An algorithm for total variation minimization and
% applications,” Journal of Mathematical Imaging and Vision, vol. 20,
% pp. 89-97, 2004.
%
%
% For further details about the TwIST algorithm, see the paper:
%
% J. Bioucas-Dias and M. Figueiredo, "A New TwIST: Two-Step
% Iterative Shrinkage/Thresholding Algorithms for Image 
% Restoration",  IEEE Transactions on Image processing, 2007.
% 
% and
% 
% J. Bioucas-Dias and M. Figueiredo, "A Monotonic Two-Step 
% Algorithm for Compressive Sensing and Other Ill-Posed 
% Inverse Problems", submitted, 2007.

%
% Paper: J. Bioucas- Dias and M.  Figueiredo, "A New TwIST: Two-Step
% Iterative Shrinkage/Thresholding Algorithms for Image Restoration", 
% Submitted to IEEE Transactions on Image processing,  2007.
%
% (available at   http://www.lx.it.pt/~bioucas/publications.html)
%
%
% 
% Authors: Jose Bioucas-Dias and Mario Figueiredo, 
% Instituto Superior Técnico, Outubro, 2007

clear all
close all



%%%%%%%%%%%%%%%%%%%%% Original Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = double(imread('cameraman.tif'));
%x = double(imread('phantom.tif'));
%x=256*phantom(256);

N=length(x);
colormap(gray);

% remove mean (not necessary)
mu=mean(x(:));
x=x-mu;

% define the blur operator
% N must be even
middle=N/2+1;
% blurr matrix
B=zeros(N);
% Uniform blur
lx=4; %blur x-size
ly=4; %blurr y-size
B((middle-ly):(middle+ly),(middle-lx):(middle+lx))=1;
%circularly center
B=fftshift(B);
%normalize
B=B/sum(sum(B));
% convolve
y = real(ifft2(fft2(B).*fft2(x)));

% set BSNR
BSNR = 40;
Py = var(y(:));
sigma= sqrt((Py/10^(BSNR/10)));

% add noise
y=y+ sigma*randn(N);

% plot figures
figure(1); colormap gray; 
imagesc(x); axis off;
title('Original image')
figure(2); colormap gray; 
title('Noisy and blurred image')
imagesc(y); axis off;


% smoothing parameter (empirical setting)
tau = 2e-2*sigma^2/0.56^2;

% extreme eigenvalues (TwIST parameter)
lam1=1e-4;    
% TwIST is not very sensitive to this parameter
% rule of thumb: lam1=1e-4 for severyly ill-conditioned% problems
%              : lam1=1e-1 for mildly  ill-conditioned% problems
%              : lam1=1    when A = Unitary matrix


% ------------  TV experiments ---------------------
K=fft2(B);
KC=conj(K);

% handle functions for TwIST
%  convolution operators
A = @(x) real(ifft2(K.*fft2(x)));
AT = @(x) real(ifft2(KC.*fft2(x)));

% denoising function;
tv_iters = 5;
Psi = @(x,th)  tvdenoise(x,2/th,tv_iters);
% TV regularizer;
Phi = @(x) TVnorm(x);
%Phi = @(x) sum(sum(sqrt(diffh(x).^2+diffv(x).^2)));


% start with the wiener filter
varx = var(y(:)); 	% approximate var of x
x0 = real(fft2(KC./(abs(KC).^2+10*sigma^2/varx).*ifft2(y)));

tolA = 1e-4;
% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'
[x_twist,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
         TwIST(y,A,tau,...
         'AT', AT, ...
         'lambda',lam1,...
         'True_x', x,...       
         'Psi', Psi, ...
         'Phi',Phi, ...
         'Monotone',1,...
         'Initialization',x0,...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 1);

figure(3); colormap gray; 
title('TwIST restored image')
imagesc(x_twist); axis off;
drawnow;


     
% -- IST (lam1=1) ---------------------------
% stop criterium:  the objective function 
% falls below obj_twist(end)
%
% IST takes too long and thus we run only 200 iterations
[x_ist,dummy,obj_ist,...
  times_ist,dummy,mse_ist]= ...
         TwIST(y,A,tau,...
         'AT', AT, ...
         'lambda',1,...
         'True_x', x,...       
         'Psi', Psi, ...
         'Phi',Phi, ...
         'Monotone',1,...
         'MaxiterA',200,...
         'Initialization',x0,...
         'StopCriterion',3,...
       	 'ToleranceA',obj_twist(end),...
         'Verbose', 1);
     
     
     
figure(4); colormap gray; 
title('IST restored image')
imagesc(x_ist); axis off;
drawnow;
     
   
figure(5)
subplot(2,1,1)
semilogy(times_twist,obj_twist,'r',...
         times_ist,obj_ist,'b','LineWidth',2)
legend('TwIST','IST')
st=sprintf('tau = %2.2e, sigma = %2.2e',tau,sigma),...
title(st)
ylabel('Obj. function')
xlabel('CPU time (sec)')

grid
subplot(2,1,2)
plot(times_twist(2:end),mse_twist(2:end),'r',...
         times_ist(2:end),mse_ist(2:end),'b','LineWidth',2)
legend('TwIST','IST')
ylabel('MSE')
xlabel('CPU time (sec)')


fprintf(1,'TwIST   CPU time - %f\n', times_twist(end));
fprintf(1,'IST     CPU time - %f\n', times_ist(end));
