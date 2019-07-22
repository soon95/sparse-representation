function varargout = dictnormalize(D)
% dictnormalize   Normalize and arrange the vectors of a dictionary (frame)
% Note that in addition to set the 2-norm of each column vector to 1
% the vectors are ordered by 'frequency', i.e. number of zero crossings.
% If only normalization 
%   use: Dn = D ./ repmat(sqrt(sum(D.^2)),[size(D,1) 1]);
%   or:  Dn = D*diag(1./sqrt(diag(D'*D)));
%   or:  Dn = D.*(ones(size(D,1),1)*(1./sqrt(sum(D.*D))));
%
% Dn = dictnormalize(D);       
% [Dn, I] = dictnormalize(D);     % Dn = Ds(:,I); where Ds = D.*(ones(K,1)*sf)
% [Dn, I, sf] = dictnormalize(D);     % sf = (+/-) 1./sqrt(sum(D.*D))
%-----------------------------------------------------------------------------------
% arguments:
%   D      the input dictionary, a NxK matrix
%   Dn     the normalized dictionary, a NxK matrix
%   I      the order of the dictionary columns, D
%   sf     the scalefactor, 1xK vector, element k is (+/-)1/norm(D(:,k))
%-----------------------------------------------------------------------------------
% ex:   N=20; K=40; L=100;
%       D = randn(N,K);  W = randn(K,L);  X = D*W;
%       [Dn, I, sf] = dictnormalize(D);
%       D1 = D(:,I).*(ones(N,1)*sf(I));     % how Dn was made
%       disp(norm(Dn-D1));   % should be 0
%       Wn = zeros(size(W)); % coefficients for normalized dictionary
%       for k=1:K; Wn(k,:) = W(I(k),:)/sf(I(k)); end;
%       Xn = Dn*Wn;
%       disp(norm(X-Xn));   % should be 0


%----------------------------------------------------------------------
% Copyright (c) 2009.  Karl Skretting.  All rights reserved.
% University of Stavanger, Signal Processing Group
% Mail:  karl.skretting@uis.no   Homepage:  http://www.ux.uis.no/~karlsk/
% 
% HISTORY:  dd.mm.yyyy
% Ver. 1.0  12.02.2009  KS: function made based on ...\FrameTools\NormalizeF 
% Ver. 1.1  13.03.2009  KS: more output arguments
%----------------------------------------------------------------------

Mfile='dictnormalize';

if (nargin ~= 1)
   error([Mfile,': wrong number of arguments, see help.']);
end

[N,K] = size(D);

% normalize each of these
sf = ones(1,K);
for k=1:K;    % now normalize Fb
   temp = D(:,k)'*D(:,k);
   if (temp > 0); sf(k) = 1/sqrt(temp); end;
   if (sum(D(:,k)) < 0); sf(k) = -sf(k); end;
   D(:,k) = D(:,k) * sf(k);
end

% set the order in decreasing 'frequency'
% t = D(2:N,:) - D(1:(N-1),:);
% t = sum(abs(t));
% [t,I] = sort(t);
% D = D(:,I);

if (nargout >= 1); varargout(1) = {D}; end;
if (nargout >= 2); varargout(2) = {I}; end;
if (nargout >= 3); varargout(3) = {sf}; end;

return