%% Script to compute PIAA lens sag profiles for KPIC
% G. Ruane
% Based on Guyon 2003, Galicher 2005, Vanderbei 2005

clear; %close all;
addpath(genpath('PIAA_lib'));

L = 10;% Distance between the PIAA lenses
lambda0 = 2.2;% Wavlength (microns)
deltaLambda = lambda0*0.2;% Bandwidth (microns)
numWavelengths = 9;% Number of wavelength samples across band
Npts = 10001;% Number of points for design
Nrays = 101; % Number of rays to send through PIAA
material = 'CaF2';


label = ['PIAAsag_',material,'_L',num2str(L),'_lam',num2str(lambda0),...
    '_delLam',num2str(deltaLambda),'_numLam',num2str(numWavelengths)];


%% Get material properties

% Define the refractive index of each lens 

% at central wavelength
n1_lam0 = getRefractiveIndex(material,lambda0);
n2_lam0 = getRefractiveIndex(material,lambda0);

% across the whole band 
lambdas = linspace(lambda0-deltaLambda/2,lambda0+deltaLambda/2,numWavelengths);
n1 = getRefractiveIndex(material,lambdas);
n2 = getRefractiveIndex(material,lambdas);

%% Get remapping function

% Annulus to gaussian remapping function
Rin1 = 0.236; 
Rout1 = 1; 
Rin2 = 0.1; 
Rout2 = 1; 
w = 1/sqrt(2);

[r1,r2] = annulus2truncGaussianRemappingPIAA(Rin1,Rout1,Rin2,Rout2,w,Npts);

%% Make the PIAA sag profiles 

% PIAA is designed for central wavelength
PIAA = makePIAAlenses(r1,r2,n1_lam0,n2_lam0,L);

%% Ray tracing 

plotPIAAraytrace_polychr;
asdf
lensThickness = 0.3;

figure;
% colorOrd = get(gca,'ColorOrder');
colorOrd = jet(numWavelengths);

for lambdaIndex = 1:numWavelengths

    % Update refractive index 
    PIAA.lens1.n = n1(lambdaIndex);
    PIAA.lens2.n = n2(lambdaIndex);

    % Trace rays
    [RAYS,PIAA] = rayTracePIAA(PIAA,Nrays,false);

    % Plot the rays 
    plot(RAYS.z,RAYS.x,'Color',colorOrd(lambdaIndex,:));hold on;

end

% Update refractive index 
PIAA.lens1.n = n1_lam0;
PIAA.lens2.n = n2_lam0;
[RAYS,PIAA] = rayTracePIAA(PIAA,Nrays,false);
plot(RAYS.z,RAYS.x,'k');
     
xLens1 = PIAA.lens1.xFull;
xLens2 = PIAA.lens2.xFull;
zLens1 = PIAA.lens1.zFull;
zLens2 = PIAA.lens2.zFull; 

% Plot the lens surfaces
plot(zLens1,xLens1,'k','LineWidth',2);
plot(zLens2,xLens2,'k','LineWidth',2);

% Plot the lens outlines (for asthetics)
plot([zLens1(end)-lensThickness zLens1(end)],[Rout1 Rout1],'k','LineWidth',2);
plot([zLens1(end)-lensThickness zLens1(end)],[-Rout1 -Rout1],'k','LineWidth',2);
plot([zLens1(end)-lensThickness zLens1(end)-lensThickness],[-Rout1 Rout1],'k','LineWidth',2);
plot([zLens2(end) zLens2(end)+lensThickness],[Rout2 Rout2],'k','LineWidth',2);
plot([zLens2(end) zLens2(end)+lensThickness],[-Rout2 -Rout2],'k','LineWidth',2);
plot([zLens2(end)+lensThickness zLens2(end)+lensThickness],[-Rout2 Rout2],'k','LineWidth',2);

hold off; 
xlabel('z/a');
ylabel('x/a');
axis equal
axis([min(RAYS.z(:)) max(RAYS.z(:)) -1.5 1.5])