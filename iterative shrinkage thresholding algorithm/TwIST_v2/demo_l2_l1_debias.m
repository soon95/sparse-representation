function demo_l2_l1_debias
% demo_l2_l1 - This demo illustrates  the TwIST  
% algorithm in  the l2-l1 optimization problem 
%
%     xe = arg min 0.5*||A x-y||^2 + tau ||x||_1
%             x
%
% where A is a generic matrix and ||.||_1 is the l1 norm.
% After obtaining the solution we implement a debias phase
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
n_spikes = 160;
% random +/- 1 signal
x = zeros(n,1);
q = randperm(n);
x(q(1:n_spikes)) = sign(randn(n_spikes,1));

% measurement matrix
disp('Building measurement matrix...');
R = randn(k,n);
%normalize R (not necessary)
% maxSingValue=svds(R,1);
% R=R/maxSingValue;
% R has been precomputed 
% load CSmatrix
% R=R;
disp('Finished creating matrix');

%TwIST handlers
% Linear operator handlers
hR = @(x) R*x;
hRt = @(x) R'*x;
% define the regularizer and the respective denoising function
% TwIST default
%Psi = @(x,th) soft(x,th);   % denoising function
%Phi = @(x) l1norm(x);       % regularizer

% noise variance
sigma=1e-2;
% observed data
y = hR(x)+sigma*randn(k,1);

% regularization parameter 
tau = 0.1*max(abs(hRt(y)));

% TwIST parameters
lambda1 = 1e-3;  

%             If min eigenvalue of A'*A == 0, or unknwon,  
%             set lam1 to a value much smaller than 1. 
%             TwIST is not very sensitive to this parameter
%             
%
%             Rule of Thumb: 
%                 lam1=1e-4 for severyly ill-conditioned problems
%                 lam1=1e-2 for mildly  ill-conditioned problems
%                 lam1=1    for A unitary direct operators


% stopping theshold
tolA = 1e-4;

% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'
[x_twist,x_debias_twist,obj_twist,...
    times_twist,debias_start_twist,mse]= ...
         TwIST(y,hR,tau, ...
         'Lambda', lambda1, ...
         'Debias',1,...
         'AT', hRt, ... 
         'Monotone',1,...
         'Initialization',0,...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 1);
    

h_fig = figure(2);
set(h_fig,'Units','characters',...
        'Position',[10 10 250 45])
subplot(2,1,1)
semilogy(times_twist(1:debias_start_twist),obj_twist(1:debias_start_twist),'r',...
    times_twist(debias_start_twist:end),obj_twist(debias_start_twist:end),'b','LineWidth',2)
legend('TwIST', 'Debias phase')
st=sprintf('\\lambda = %2.1e, sigma = %2.1e',tau,sigma);
title(st)
xlabel('CPU time (sec)')
ylabel('Obj. function')
grid
subplot(2,1,2)
plot(1:n,x_debias_twist,'b', 1:n, x+2.5,'k','LineWidth',2)
axis([1 n -1, 5]);
legend('TwIST','Original','x')
st=sprintf('TwIST MSE = %2.1e', sum((x_debias_twist-x).^2)/prod(size(x)));
title(st)

fprintf('TwIST  MSE = %2.1e\n', sum((x_debias_twist-x).^2)/prod(size(x)));
fprintf(1,'TwIST   CPU time - %f\n', times_twist(end));

