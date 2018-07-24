function [ out ] = lpgauss(H,W,SIGMA)

%   Creates a 2D Gaussian filter for a Fourier space image
%   W is the number of columns of the source image and H is the number of
%   rows. SIGMA is the standard deviation of the Gaussian.

H = single(H);
W = single(W);
kcx = (SIGMA);
kcy = ((H/W)*SIGMA);
[x,y] = meshgrid(-floor(W/2):floor((W-1)/2), -floor(H/2):floor((H-1)/2));
out = ifftshift(exp(-(x.^2/(kcx^2)+y.^2/(kcy^2))));

end