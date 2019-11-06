function y = masked_FFT_t(x,mask)
%
% This is the transpose operator of the
% partial FFT transform, implemented in 
% masked_FFT.m.  See that file for details.
%
[n1,n2] = size(mask);
if n1~=n2
   disp('only square images please!')
   return
end
n = n1;
if floor(log2(n))~=log2(n)
   disp('only power of 2 sizes, please!')
   return
end
ii=find(abs(mask)>0);
gg = zeros(n,n); 
gg(ii)=x./mask(ii);
y = real(ifft2(ifftshift(gg)))*n;



