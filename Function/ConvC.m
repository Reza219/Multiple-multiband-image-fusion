function [A] = ConvC(X,FK,nr)
% defines a circular convolution (the same for all bands)
% FK is the fft2 of a one-band filter 

[Ns,Np]=size(X);
nc = Np/nr;
A = reshape(real(ifft2(fft2(reshape(X',[nr,nc,Ns])).*repmat(FK,[1,1,Ns]))),nr*nc,Ns)';