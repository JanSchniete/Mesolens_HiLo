function [ out ] = bpgauss( H,W,SIGMA )

%   Input height(H), width(W) and standard deviation (SIGMA) to construct
%   the bandpass filter

H = single(H);
W = single(W);
kcx = (SIGMA);
kcy = ((H/W)*SIGMA);
[x,y] = meshgrid(-floor(W/2):floor((W-1)/2), -floor(H/2):floor((H-1)/2));
out = ifftshift(exp(-(x.^2/(2*kcx^2)+y.^2/(2*kcy^2)))-exp(-(x.^2/(kcx^2)+y.^2/(kcy^2))));

end

