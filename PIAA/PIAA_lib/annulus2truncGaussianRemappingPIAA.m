function [r1,r2] = annulus2truncGaussianRemappingPIAA(Rin1,Rout1,Rin2,Rout2,a,Npts)
%[r1,r2] = annulus2truncGaussianRemappingPIAA(Rin1,Rout1,Rin2,Rout2,a,Npts)
%   Generates pairs of radii that maps points from an evenly illuminated
%   annular pupil to a 'truncated' Gaussian with std dev of sigma 
%   
%   Inputs:
%       Rin1  - Inner radius of annulus (zero for circular pupil)
%       Rout1 - Outer radius of annulus  
%       Rin2  - Inner radius of the truncated Gaussian 
%       Rout2 - Outer radius of the truncated Gaussian   
%       a - Output Gaussian width
%       Npts - Number of points (linearly spaced betweern Rin and Rout)
%
%   Outputs:
%       r1,r2 - Pairs of radii that a single ray intersects each lens.

    r1 = linspace(Rin1, Rout1, Npts);% radii of input rays 
    
    % Conservation of encircled energy (see mathematica notebook)
    numer = exp((Rin2^2 + Rout2^2)/a^2)*(Rout1^2-Rin1^2);
    
    denom = exp(Rin2^2/a^2).*(r1 - Rin1).*(r1 + Rin1) + ...
                exp(Rout2^2/a^2)*(Rout1^2 - r1.^2);
            
    r2 = a*sqrt(log(numer./denom));
    
end

