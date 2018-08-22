%% Script to compute PIAA lens sag profiles for KPIC
% G. Ruane
% Based on Guyon 2003, Galicher 2005, Vanderbei 2005

clear; %close all;
addpath(genpath('PIAA_lib'));

L = 10;% Distance between the PIAA lenses
lambda = 2.2;% Wavlength (microns)
Npts = 10001;% Number of rays for design 
material = 'CaF2';

label = ['PIAAsag_',material,'_L',num2str(L),'_lam',num2str(lambda)];

%% Get material properties

n1 = getRefractiveIndex(material,lambda);
n2 = getRefractiveIndex(material,lambda);

%% Get remapping function

% Annulus to gaussian remapping function
Rin1 = 0.236; 
Rout1 = 1; 
Rin2 = 0.1; 
Rout2 = 1; 
w = 1/sqrt(2);

[r1,r2] = annulus2truncGaussianRemappingPIAA(Rin1,Rout1,Rin2,Rout2,w,Npts);

%% Make the PIAA sag profiles 

PIAA = makePIAAlenses(r1,r2,n1,n2,L);

%% Ray tracing 
Nrays = 21;

[RAYS,PIAA] = rayTracePIAA(PIAA,Nrays,false);


%% Make Plot
plotPIAAraytrace;
asdf 

xAxisLimit = 1.5*max([Rout1 Rout2]);

lensThickness = 0.3;

xLens1 = PIAA.lens1.xFull;
xLens2 = PIAA.lens2.xFull;
zLens1 = PIAA.lens1.zFull;
zLens2 = PIAA.lens2.zFull; 

figure;
    colorOrd = get(gca,'ColorOrder');
    
    % Plot the lens surfaces
    plot(zLens1,xLens1,'k','LineWidth',2);hold on;
    plot(zLens2,xLens2,'k','LineWidth',2);
    
    % Plot the lens outlines (for asthetics)
    plot([zLens1(end)-lensThickness zLens1(end)],[Rout1 Rout1],'k','LineWidth',2);
    plot([zLens1(end)-lensThickness zLens1(end)],[-Rout1 -Rout1],'k','LineWidth',2);
    plot([zLens1(end)-lensThickness zLens1(end)-lensThickness],[-Rout1 Rout1],'k','LineWidth',2);
    plot([zLens2(end) zLens2(end)+lensThickness],[Rout2 Rout2],'k','LineWidth',2);
    plot([zLens2(end) zLens2(end)+lensThickness],[-Rout2 -Rout2],'k','LineWidth',2);
    plot([zLens2(end)+lensThickness zLens2(end)+lensThickness],[-Rout2 Rout2],'k','LineWidth',2);
    
    % Plot the rays 
    plot(RAYS.z,RAYS.x,'Color',colorOrd(1,:));
    
    xlabel('z/a');
    ylabel('x/a');
    axis equal
    axis([min(RAYS.z(:)) max(RAYS.z(:)) -xAxisLimit xAxisLimit])