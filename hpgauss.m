function [ out ] = hpgauss(H,W,SIGMA)

%   Creates a 2D Gaussian filter for a Fourier space image of height H and
%   width W. SIGMA is the standard deviation of the Gaussian.

out=1-lpgauss(H,W,SIGMA);

end

