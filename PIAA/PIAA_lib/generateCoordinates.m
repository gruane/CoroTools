function coords = generateCoordinates( N )
%coords = generateCoordinates( N )
%   Generates sample centered coordinate system (both cartesian and polar)

    % Create coordinate system 
    [X,Y] = meshgrid(-N/2:N/2-1);
    [THETA,RHO] = cart2pol(X,Y);
    xvals = X(1,:);yvals = Y(:,1);
    
    coords.X = X;
    coords.Y = Y;
    coords.THETA = THETA;
    coords.RHO = RHO;
	coords.xvals = xvals;
    coords.yvals = yvals;
end

