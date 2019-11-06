function demo_Piecewise_cubic_polynomial
%
% In this problem, the TwIST algorithm introduced in
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
% solves the l2-l1 optimization problem 
%
%     xe = arg min 0.5*||A x-y||^2 + tau ||x||_1,
%             x
% where ||.||_2 and ||.||_1 are the l2 and l1 norms, respectively,  and 
% the linear operator A and the data vector y correspond  to  the 
% piecewise cubic polynomial sparse  reconstruction problem stated in 
%
% E. J. Candés and J. Romberg. Practical signal recovery from random 
% projections. In Wavelet Applications in Signal and Image Processing XI, 
% Proc. SPIE Conf. 5914., 2004.
%
% -------------------------------------------------------------------------
%  SPARCO toolbox has been used to build the problem environment. For
%  details,  we refer to 
%
%  E. van den Berg and M. P. Friedlander and G. Hennenfent and F. Herrmann
%  and R. Saab and O. Yilmaz, "Sparco: testing framework for sparse
%  reconstruction", University of British Columbia, Dept. Computer Science,
%  Vancouver, 2007.
%
%  http://www.cs.ubc.ca/labs/scl/sparco/
%
% -------------------------------------------------------------------------
%
% Please check for the latest version of the code and papers at
%
% www.lx.it.pt/~bioucas/TwIST
%
%
% Authors: Jose Bioucas-Dias and Mario Figueiredo, 
% Instituto Superior Técnico, October, 2007
%
 
% Generate environment
  P = generateProblem(6);
  
  A  = @(x) P.A(x,1);      % The direct operator
  AT = @(y) P.A(y,2);      % and its transpose.

  b  = P.b;                % Data vector.
  
  m = P.sizeA(1);          % m is the no. of rows.
  n = P.sizeA(2);          % n is the no. of columns.
  
% Solve an L1 recovery problem:
% minimize  1/2|| Ax - b ||_2^2  +  tau ||x||_1.
  
   tau = 100;    % regularization parameter
%   [x,dummy,obj,time] = GPSR_BB(b, A, tau, ...
%      'AT', AT, ...
%      'Continuation',1, ...
%      'Monotone', 1, ...
%      'StopCriterion',1,...
%      'ToleranceA',1e-5);
 
  [x,dummy,obj,time] = TwIST(b, A, tau, ...
     'AT', AT, ...
     'Monotone', 1, ...
     'StopCriterion',1,...
     'ToleranceA',1e-5);

  
% The solution x is the reconstructed signal in the sparsity basis. 
  h_fig=figure(1);
  set(h_fig,'Units','characters',...
        'Position',[0 35 100 30]);
  plot(x,'LineWidth',2);
  title('Coefficients of the reconstructed signal')
  xlabel('Coefficient')
  
% Use the function handle P.reconstruct to use the coefficients in
% x to reconstruct the original signal.
  y     = P.reconstruct(x);    % Use x to reconstruct the signal.
  yorig = P.signal;            % P.Signal is the original signal.

  h_fig=figure(2);
  set(h_fig,'Units','characters',...
        'Position',[100 35 100 30]);

  plot(1:length(y), y    ,'b', ...
       1:length(y), yorig,'r','LineWidth',2);
  legend('Reconstructed signal','Original signal');
  title('Reconstructed and original signals');

% plot the evolution of the objective function
  h_fig=figure(3);
  set(h_fig,'Units','characters',...
        'Position',[0 0 100 30]);

  semilogy(time,obj,'LineWidth',2)
  title('Evolution of the objective function');  
  xlabel('CPU time (sec)');
  grid on
 

  
end % function demo_Piecewise_cubic_polynomial
