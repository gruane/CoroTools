function [RAYS,PIAA] = rayTracePIAA(PIAA,Nrays,showPlot)
%[RAYS,PIAA] = rayTracePIAA(PIAA,Nrays,showPlot)
%   Launches Nrays rays through the PIAA optics and returns the
%   intersection points as two arrays RAYS.z and RAYS.x. 
%   
%   Inputs:
%       PIAA - The PIAA structure from makePIAAlenses function.
%       Nrays - The number of rays for the ray trace. 
%       showPlot - True to show a plot 
%
%   Outputs:
%       RAYS - Structure with (x,z) points of the ray intersections.
%       PIAA - Updated PIAA structure with PIAA reflected over the z-axis.

    L = PIAA.L; % Distance between the PIAA lenses
    n1 = PIAA.lens1.n;% refractive index of lens 1
    n2 = PIAA.lens2.n;% refractive index of lens 2 
    
    a1 = max(PIAA.lens1.r);% Radius of lens 1 (for plotting)
    a2 = max(PIAA.lens2.r);% Radius of lens 2 (for plotting)

    % Allocation of the output arrays of ray points  
    RAYS.z = zeros(4,Nrays);
    RAYS.x = zeros(4,Nrays);

    % Mirror the lens profiles over the z-axis 
    xLens1 = [-1*fliplr(PIAA.lens1.r(2:end)) PIAA.lens1.r];
    xLens2 = [-1*fliplr(PIAA.lens2.r(2:end)) PIAA.lens2.r];
    zLens1 = [fliplr(PIAA.lens1.z(2:end)) PIAA.lens1.z];
    zLens2 = [fliplr(PIAA.lens2.z(2:end)) PIAA.lens2.z];

    % Returns the full PIAA surface used for ray tracing 
    PIAA.lens1.xFull = xLens1;
    PIAA.lens2.xFull = xLens2;
    PIAA.lens1.zFull = zLens1;
    PIAA.lens2.zFull = zLens2;  
    
    % Define launch and termination position along z axis for rays 
    zRayLaunch = -L;
    zRayFinish = 2*L;
    
    % Lens thickness for plot
    lensThickness = 0.05*L;

    % Get the local slopes of the lenses 
    zpLens1      = gradient(zLens1)./gradient(xLens1);
    zpLens2      = gradient(zLens2)./gradient(xLens2);

    x0s = linspace(-1,1,Nrays);% Ray launch points in transverse direction
    zIntLens1 = interp1(xLens1,zLens1,x0s);% z-position where ray intersects lens 1
    zpIntLens1 = interp1(xLens1,zpLens1,x0s);% the slope of the lens where ray intersects lens 1
    alpha = atan(1./zpIntLens1);% Surface angle (radians)

    theta1Lens1 = pi/2-alpha;% The angle of incidence of the rays at lens 1 (radians)
    theta1Lens1(theta1Lens1>pi/2) = theta1Lens1(theta1Lens1>pi/2) - pi;
    theta2Lens1 = asin(n1*sin(theta1Lens1));% Refraction angle of the rays after lens 1 (radians)
    raySlopes = tan(theta2Lens1+alpha-pi/2);% slopes of the rays after lens 1 

    % Initial (z,x) position of the rays 
    RAYS.z(1,:) = zRayLaunch;
    RAYS.x(1,:) = x0s;
    
    % (z,x) position of the rays at lens 1 surface 
    RAYS.z(2,:) = zIntLens1;
    RAYS.x(2,:) = x0s;

    % Make a plot if requested 
    if(showPlot)
        showNorms = false;% Shows the surface normals 

        figure(1);
        colorOrd = get(gca,'ColorOrder');
        plot(zLens1,xLens1,'k','LineWidth',2);hold on;
        plot(zLens2,xLens2,'k','LineWidth',2);
        plot([zLens1(end)-lensThickness zLens1(end)],[a1 a1],'k','LineWidth',2);
        plot([zLens1(end)-lensThickness zLens1(end)],[-a1 -a1],'k','LineWidth',2);
        plot([zLens1(end)-lensThickness zLens1(end)-lensThickness],[-a1 a1],'k','LineWidth',2);
        plot([zLens2(end) zLens2(end)+lensThickness],[a2 a2],'k','LineWidth',2);
        plot([zLens2(end) zLens2(end)+lensThickness],[-a2 -a2],'k','LineWidth',2);
        plot([zLens2(end)+lensThickness zLens2(end)+lensThickness],[-a2 a2],'k','LineWidth',2);
    end
    
    % Loop over rays 
    for rayNum = 1:Nrays

        % Get the ray from lens 1 to lens 2
        zRayTmp = linspace(zIntLens1(rayNum),zRayFinish,1000);
        ray1to2 = raySlopes(rayNum)*(zRayTmp-zIntLens1(rayNum)) + x0s(rayNum);

        % Find the intersection between the ray and lens 2
        [zIntLens2,xIntLens2] = intersections(zRayTmp,ray1to2,zLens2,xLens2,false);

        % Get the slope of lens 2 at the intersection
        zpIntLens2 = interp1(xLens2,zpLens2,xIntLens2);% the slope of the lens where ray intersects lens 2
        
        % Apply Snell's law at the second surface 
        beta = atan(1./zpIntLens2);% Surface angle (radians)
        theta1Lens2 = atan(raySlopes(rayNum))-beta+pi/2;% The angle of incidence of the rays at lens 2 (radians)
        if(theta1Lens2>pi/2);theta1Lens2 = theta1Lens2 - pi; end
        theta2Lens2 = asin(1/n2*sin(theta1Lens2));
        raySlopeOut = tan(theta2Lens2+beta-pi/2);% slopes of the rays after lens 2 

        % Get termination point 
        finalXpoint = raySlopeOut*(zRayFinish-zIntLens2) + xIntLens2;
        
        % Output ray points along surface 2
        RAYS.z(3,rayNum) = zIntLens2;
        RAYS.x(3,rayNum) = xIntLens2;
        
        % Output ray points at termination of rays 
        RAYS.z(4,rayNum) = zRayFinish; % Last z point 
        RAYS.x(4,rayNum) = finalXpoint;
        
        % Make a plot if requested 
        if(showPlot)
            
            plot([zRayLaunch zIntLens1(rayNum)],[x0s(rayNum) x0s(rayNum)],'Color',colorOrd(1,:));
            
            if(showNorms)
                % Plot the normal to lens 1 at the intersection
                normLine = -1*zpIntLens1(rayNum)*(zLens1-zIntLens1(rayNum)) + x0s(rayNum);
                plot(zLens1,normLine,'r');
            end

            % Plot the ray from lens 1 to lens 2
            plot([zIntLens1(rayNum) zIntLens2],interp1(zRayTmp,ray1to2,[zIntLens1(rayNum) zIntLens2]),'Color',colorOrd(1,:));

            if(showNorms)
                % Plot the normal to lens 2 at the intersection
                normLine = -1*zpIntLens2*(zLens2-zIntLens2) + xIntLens2;
                plot(zLens2,normLine,'r');
            end

            plot([zIntLens2 zRayFinish],[xIntLens2 finalXpoint],'Color',colorOrd(1,:));
        end
        
    end
    
    % Finish off the plot 
    if(showPlot)
        hold off;
        xlabel('z/a');
        ylabel('x/a');
        axis equal
        axis([zRayLaunch zRayFinish -1.5 1.5])
    end

end

