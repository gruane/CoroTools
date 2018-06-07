function w = generateTukeyWindow( Nwindow, RHO, alpha )
%generateTukeyWindow Generates a 2D Tukey window
%   Inputs: 'Nwindow' - the size of the window in samples
%           'RHO' - a 2D array of radial coordinates 
%           'alpha' - the squareness factor 
%   Outputs: 'w' - the 2D Tukey window 

    Nlut = round(10*Nwindow); % Number of samples for the lookup table 
    p = linspace(-Nwindow/2,Nwindow/2,Nlut);% Array of coordinates in samples
    lut = tukeywin(Nlut,alpha);% Returns a 1D Tukey window

    w = interp1(p,lut,RHO,'linear',0);% Generates the 2D window
    
end

