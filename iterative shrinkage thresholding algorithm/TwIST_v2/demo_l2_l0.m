function demo_l2_l0
% demo_l2_l0 - This demo illustrates  the TwIST algorithm in  
% the l2-l0 optimization problem 
%
%     xe = arg min 0.5*||A x-y||^2 + tau ||x||_0
%             x
%  
% which appears, for example, in compressed senseng (CS) 
% applications
%
% For further details about the TwIST algorithm, see the paper:
%
% J. Bioucas-Dias and M. Figueiredo, "A New TwIST: Two-Step
% Iterative Shrinkage/Thresholding Algorithms for Image 
% Restoration",  IEEE Transactions on Image processing, 2007.
% 
%% 
% Please check for the latest version of the code and papers at
% www.lx.it.pt/~bioucas/TwIST
%
% Authors: Jose Bioucas-Dias and Mario Figueiredo, 
% Instituto Superior Técnico, October, 2007
%
%

% signal length 
n = 4096;
% observation length 
k = 1024;

% number of spikes
n_spikes = 100;
% random +/- 1 signal
x = zeros(n,1);
q = randperm(n);
x(q(1:n_spikes)) = sign(randn(n_spikes,1));

% measurement matrix
disp('Building measurement matrix...');
R = randn(k,n);
%normalize R
%  maxSingValue=svds(R,1);
%  R=R/maxSingValue;
% R has been precomputed 
% load CSmatrix R
disp('Finished creating matrix');

%TwIST handles
% Linear operator handles 
hR = @(x) R*x;
hRt = @(x) R'*x;
% define the regularizer and the respective denoising function
Psi = @(x,th) hard(x,th);   % denoising function
Phi = @(x) l0norm(x);       % regularizer

% noise variance
sigma=  1e-2;
% observed data
y = hR(x)+sigma*randn(k,1);

%  regularization parameter 
tau = 20;

% stopping theshold
tolA = 1e-3;

% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below tolA
[x_twist,x_debias_twist,obj_twist,...
    times_twist,debias_start_twist,mse]= ...
         TwIST(y,hR,tau,...
         'Psi',Psi,...
         'Phi',Phi,...
         'AT', hRt, ... 
         'Initialization',0,...
         'Monotone', 1, ...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 1);



figure(2)
subplot(1,2,1)
semilogy(times_twist,obj_twist,'b','LineWidth',2)
legend('TwIST')
st=sprintf('\\lambda = %2.2e, sigma = %2.2e',tau,sigma),...
title(st)
ylabel('Obj. function')
xlabel('CPU time (sec)')
grid
subplot(1,2,2)
plot(1:n,x_twist,'b', 1:n, x+2.5,'k','LineWidth',2)
axis([1 n -1.5 4 ]);
legend('TwIST', 'x')
st=sprintf('TwIST MSE = %2.2e', sum((x_twist-x).^2)/prod(size(x)));

title(st)

fprintf(1,'TwIST    MSE = %e\n', sum((x_twist-x).^2)/prod(size(x)));
fprintf(1,'TwIST   CPU time - %f\n', times_twist(end));


