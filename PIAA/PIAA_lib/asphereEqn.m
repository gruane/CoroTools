function sag = asphereEqn(x,r)
%sag = asphereEqn(x,r)
%   Asphere equation 
%   See: https://en.wikipedia.org/wiki/Aspheric_lens
    
    % x = [R, kappa, z0, alpha4, alpha6, ...]
    R = x(1); % Radius of curvature 
    kappa = x(2);% Conic constant 
    z0 = x(3);% Offset in z direction
    
    % Conic term 
    sag = z0 + r.^2./(R*(1+sqrt(1-(r/R).^2.*(1+kappa))));
    
    % Aspheric coefficients 
    for alphaIndex = 1:numel(x)-3
        expon = 2+2*alphaIndex;
        sag = sag + x(3+alphaIndex)*r.^expon;
    end
    
end

