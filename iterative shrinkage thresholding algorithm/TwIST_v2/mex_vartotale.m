function [img_estimated]=mex_vartotale(img_noisy,lambda,varargin)
% function   [img_estimated]=mex_vartotale(img_noisy,lambda,varargin)
%
%  
%        
% ========================== INPUT PARAMETERS (required) =================
% Parameter     Values
% name          and description
% ========================================================================
% img_noisy		(double) Noisy image of size ny. 
% lambda		Regularization parameter (which is multiplied by the
%               TV penalty).
%
% ======================== OPTIONAL INPUT PARAMETERS ====================
% Parameter     Values
% name          and description
% =======================================================================
% itmax         The maximum number of iterations.
%               Default: 300 
% epsilon       (double) Maximum error tolerance stopping criterium.
%               Default: 0.0005
% tau           (double) Algorithm tau parameter.
%               Default: 0.249
%
% ====================== Output parameters ===============================
% img_estimated     Estimated image
%
%