function polartrans = polarTransform( t , rstep, thetastep_deg )
%polarTransform Performs a 2D polar transform of an image with nearest
%neighbors interpolation.
%   Inputs: 't' - 2D image to transform
%           'rstep' - Radial sample spacing in transform
%           'thetastep_deg' - Azimuthal sample spacing in degrees 
%   Outputs: 'polartrans' - The polar transform of t

    [N,~] = size(t);
    t = double(t);

    center = N/2 + 1;

    thetastep = thetastep_deg*pi/180;

    polartrans = zeros(N/2/rstep,360/thetastep_deg);

    rs = 0:rstep:(N/2-rstep);
    qs = 0:thetastep:(2*pi-thetastep);

    [Rs,Qs] = meshgrid(rs,qs);

    Ylocs = center + Rs.*sin(Qs);
    Xlocs = center + Rs.*cos(Qs);

    Xlocs = round(Xlocs);
    Ylocs = round(Ylocs);
    Xlocs(Xlocs == N + 1) = N;
    Ylocs(Ylocs == N + 1) = N;

    for i = 1:numel(rs)
        for j = 1:numel(qs)
            polartrans(i,j) = t(Ylocs(j,i),Xlocs(j,i));
        end
    end

end

