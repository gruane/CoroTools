function [r1,r2] = gaussianRemappingPIAA(a1,a2,sigma,Npts)
%[r1,r2] = gaussianRemappingPIAA(a1,a2,sigma)
%   Generates pairs of radii that maps points from an evenly illuminated
%   circular pupil to a Gaussian with std dev of sigma 
%   
%   Inputs:
%       a1 - Radius of input beam 
%       a1 - Radius of output beam 
%       sigma - Beam waist size (std dev at the output lens)
%
%   Outputs:
%       r1,r2 - Pairs of radii that a single ray intersects each lens.


    r1 = linspace(0,a1,Npts);% radii of input rays 
    c = a1/sqrt(1-exp(-(a2/sigma)^2));% Norm factor
    r2 = imag(sigma*sqrt(log(1-(r1/c).^2)));% radii of output rays 
    
end

