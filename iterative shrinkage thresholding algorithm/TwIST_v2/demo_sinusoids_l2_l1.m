function demo_complex_l2_l1
% demo_complex_l2_l1 - This demo illustrates  the TwIST algorithm 
% applied to the recovering of  real sinusoids in noise. 
%
% The observation model is assumed to be 
%
%  y[k]= sum_i {a_i*cos(2*pi*f_i*k+phi_i)+n_i
%
% where n_i denotes noise, which can be represented in matrix notation as 
%   
%       y = Ax+w
%
%  where A=[Ac Ac*] and  the  columns of Ac  are a dictionary of 
% sampled complex sinusoids containing the 
% frequencies in the above summation.  With this representation,
% any component of vector x selecting a column of Ac must have a 
% a conjugate conterpart selecting the correspondent component of 
% Ac*. We dropp, however, this constraint and  and solve the l2-l1 
% optimization  problem
%
%     xe = arg min 0.5*||A x-y||^2 + tau ||x||_1.
%             x
%  
% If the vector x is sparse, the solution should  reflect the above 
% constraint.
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
%% 
% Please check for the latest version of the code and papers at
%
% www.lx.it.pt/~bioucas/TwIST
%
% Authors: Jose Bioucas-Dias and Mario Figueiredo, 
% Instituto Superior Técnico, October, 2007
%

I=sqrt(-1);

% signal length 
n = 1024;
% observation length 
k = 256;


% number of spikes
n_spikes = 10;
% spikes of amplitude 1 and uniform phase in [-pi,pi[
x = zeros(n,1);
q = randperm(n);
x(q(1:n_spikes)) = sign(randn(n_spikes,1)).* ...
                         exp(I*2*pi*rand(n_spikes,1));
                     
x=[x;conj(x)];

% measurement matrix
disp('Building measurement matrix...');
% define frequencies
W=(0:n-1)*2*pi/2/n;

% define time observations
T=(0:k-1)';
R=exp(I*T*W);

% make [A A*]
R = [R conj(R)];

%normalize R
% maxSingValue=svds(R,1);
% R=R/maxSingValue;
% % R has been precomputed 
%load CSmatrix
disp('Finished creating matrix');

%TwIST handlers
% Linear operator handlers
hR = @(x) R*x;
hRt = @(x) R'*x;
% define the regularizer and the respective denoising function
% TwIST default
Psi = @(x,th) soft(x,th);   % denoising function
Phi = @(x)    sum(abs(x(:)));     % regularizer

% noise variance
sigma=1e-3;
% observed data
y = real(hR(x))+sigma*randn(k,1);

% regularization parameter 
tau=0.5*max(abs(hRt(y))); %  strong regularization


% stopping theshold
tolA = 1e-5;

% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'
[x_twist,x_debias_twist,obj_twist,...
    times_twist,debias_start_twist,mse]= ...
        TwIST(y,hR,tau,...
         'Debias', 1 , ...
         'AT', hRt, ... 
         'Initialization',randn(size(x)),...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 1);


%%%%%%%%%%%%%%%%%%% plot results %%%%%%%%%%%%%%%%%%%%%%%%%%

h_fig=figure(2);
  set(h_fig,'Units','characters',...
        'Position',[0 10 150 45]);
subplot(2,2,1)
semilogy(times_twist(1:debias_start_twist-1),obj_twist(1:debias_start_twist-1), 'r',...
        times_twist(debias_start_twist-1:end),obj_twist(debias_start_twist-1:end),'b','LineWidth',2)
axis([times_twist(1) times_twist(end) obj_twist(debias_start_twist-1) obj_twist(1)]);
legend('TwIST','Debias phase');
st=sprintf('\\tau = %2.1e, \\sigma = %2.1e',tau,sigma);
title(st)
xlabel('CPU time (sec)')
ylabel('Obj. function')
grid
subplot(2,2,2)
plot(1:2*n,real(x), 'b', 1:2*n, imag(x)+2,'r','LineWidth',2)
axis([1 2*n -1, 4]);
legend('real(x)', 'imag(x)')
st=sprintf('Original data (number of spikes = %i)',n_spikes );
title(st)


figure(2)
subplot(2,2,3)
plot(1:2*n,real(x_debias_twist), 'b',1:2*n,real(x)+2, 'r','LineWidth',2)
axis([1 2*n -1, 4]);
legend('real(x\_twist)', 'real(x\_orig)')
st=sprintf('TwIST --- MSE (real) = %2.1e', sum(abs(real(x_debias_twist-x)).^2)/prod(size(x)));
title(st)



figure(2)
subplot(2,2,4)
plot(1:2*n,imag(x_debias_twist), 'b',1:2*n,imag(x)+2.5, 'r','LineWidth',2)
axis([1 2*n -1, 4]);
legend('imag(x\_twist)', 'imag(x\_orig)')
st=sprintf('TwIST --- MSE (imag) = %2.1e', sum(abs(imag(x_debias_twist-x)).^2)/prod(size(x)));
title(st)

fprintf('TwIST  MSE = %2.1e\n', sum((abs(x_debias_twist-x)).^2)/prod(size(x)));
fprintf(1,'TwIST   CPU time - %f\n', times_twist(end));


