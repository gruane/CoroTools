function PIAAlenses = makePIAAlenses(r1,r2,n1,n2,L)
%PIAAlenses = makePIAAlenses(r1,r2,n1,n2,L)
%   Computes sag profiles for the PIAA lens to map rays from radial distance
%   r1 to r2. The lenses have refractive indices n1 and n2. L is the
%   distance between the lenses
%   
%   Inputs:
%       r1, r2 - Pairs of radii that a single ray intersects each lens.
%                Also known as the remapping function.
%       n1, n2 - Refractive indices of the lenses, respectively. 
%       L - The distance between the two lenses (same units as r1 and r2).
%
%   Outputs:
%       PIAAlenses - Structure with (r,z) points of the lens surfaces.
%       PIAAlenses.lens1.z - The sag profile of lens 1 (same units as r1 and r2).
%       PIAAlenses.lens2.z - Same as previous, but for lens 2. 

    PIAAlenses.lens1.r = r1;
    PIAAlenses.lens1.n = n1;
    PIAAlenses.lens2.r = r2;
    PIAAlenses.lens2.n = n2;
    
    PIAAlenses.L = L;
    
    A = r1 - r2;
    
    z1 = 0;
    z2 = L;

    for rayNum = 1:numel(r1)-1

        B = z2(rayNum) - z1(rayNum);

        slopeOfSurfNorm1 = A(rayNum)/(n1*sqrt(A(rayNum)^2+B^2)-B); % Galicher 2005
        slopeOfSurfNorm2 = A(rayNum)/(n2*sqrt(A(rayNum)^2+B^2)-B); % Galicher 2005

        z1(rayNum+1) = -1*slopeOfSurfNorm1*(r1(rayNum+1)-r1(rayNum)) + z1(rayNum);
        z2(rayNum+1) = -1*slopeOfSurfNorm2*(r2(rayNum+1)-r2(rayNum)) + z2(rayNum);

    end

    PIAAlenses.lens1.z = z1;
    PIAAlenses.lens2.z = z2;
    
end

