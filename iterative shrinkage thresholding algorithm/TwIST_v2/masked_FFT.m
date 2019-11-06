function y = masked_FFT(x,mask)
%
% Computes a 2D partial FFT transform,
% only for those frequencies for which 
% mask is not zero. The result is returned
% in an array of size (k,1) where k is
% the number of non-zeros in mask.
% The transpose of this operator is 
% available in the function masked_FFT_t.m
%
% Copyright, 2006, Mario Figueiredo.
% 
%
[n1 n2] = size(x);
if n1~=n2
   disp('only square images, please!')
   return
end
n = n1;
if floor(log2(n))~=log2(n)
   disp('only power of 2 sizes, please!')
   return
end
if sum(sum(size(mask)~=size(x)))~=0
   disp('Mask must have the same size as x');
   return
end
ii = find(abs(mask)>0);
Rf = fftshift(fft2(x)).*mask/n;
y = Rf(ii);


