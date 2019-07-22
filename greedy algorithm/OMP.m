function [new_W,new_Gamma,err] = OMP(signal,dico,maxIter)
% OMP Naive implementation of Orthogonal Matching Pursuit with fixed
% dictionaries
% Inputs : - signal  the n by 1 data
%          - dico    the dictionary n by m in which to decompose the signal
%          - maxIter the maximum number of iterations

% Outputs : - W        weights of the selected atoms
%           - Gamma    set of selected atoms indexes
%           - error    vector of achieved errors in dB


% Manuel Moussallam - manuel.moussallam@telecom-paristech.fr

%-------------------- Initialization --------------------------------------
R = signal;                             % The residual R^0 f
approx = zeros(length(signal),1);       % Current approximation f_o
W = zeros(maxIter,1);                   % Weights of selected atoms
Gamma = zeros(maxIter,1);               % Set of indexes of selected atoms
error = zeros(maxIter,1);               % vector of error values  
err=[R'];
%--------------------------------------------------------------------------

% As reference, compute energy of original signal
sigEnergy = 10*log10(norm(signal));

% Main MP Loop: only max iteration number stopping criterion implemented
% here
for it=1:maxIter
   
    %---- Step 1 : Selection Step -----------%
    projs = R*dico;                   % projecting residual onto dictionary
    [w Gamma(it)] = max(abs(projs));  % selecting atom as arg max 
    
    %---- Step 2 : Update Step -------------%
    
    % Orthogonal projection on the subspace spanned by all 
    % previously selected atoms
    W(1:it) = pinv(dico(:,Gamma(1:it)))*signal'; %
    
    % update approximant and residual
    approx = dico(:,Gamma(1:it))*W(1:it);  
    pre_R=R;
    R = signal - approx';
    
    err=[err R'];
    
    if norm(R-pre_R)/length(R)<5e-4
        break;
    end
    
end

for i=1:maxIter
    if Gamma(i)~=0
        new_Gamma(i)=Gamma(i);
    else
        break;
    end
end

for i=1:maxIter
    if W(i)~=0
        new_W(i)=W(i);
    else
        break;
    end
end

end
% ------ end of main Loop -----------------%

% residual = signal;
% approx = zeros(length(residual),1);
% coeffs = zeros(maxIter,1);
% indexes =zeros(maxIter,1);
% SRR = zeros(maxIter,1);
% sigEnergy = 10*log10(norm(signal));
% for it=1:maxIter
%    
%     % Selection Step
%     projs = residual*dico;
%     [coeff indexes(it)] = max(abs(projs));
%     
%     coeffs(it) = projs(indexes(it));
%     
%     % projection Step
%     projScores = pinv(dico(:,indexes(1:it)))*signal';
%     approx = dico(:,indexes(1:it))*projScores;
%     residual = signal - approx';
%     
%     
%     
%     % srr computation
%     SRR(it) = srr(approx,residual)- sigEnergy;
% %     SRR(it) = 10*log10(norm(approx)/norm(residual));
%     
%     %escape loop
%     if SRR(it) < -120
%         disp(['Exact reconstruction found in ' num2str(it) ' iterations']);
%         SRR(it+1:end) = SRR(it);
%         break; 
%     end
% 
% end

% end

