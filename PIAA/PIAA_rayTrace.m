%% Script to compute PIAA lens sag profiles for KPIC
% G. Ruane
% Based on Guyon 2003, Galicher 2005, Vanderbei 2005

clear; close all;
addpath(genpath('PIAA_lib'));

L = 5;% Distance between the PIAA lenses
lambda = 2.2;% Wavlength (microns)
Npts = 10001;% Number of points for design
Nrays = 21; % Number of rays to send through PIAA
material = 'CaF2';

label = ['PIAAsag_',material,'_L',num2str(L),'_lam',num2str(lambda),'_Npts',num2str(Npts)];

%% Get material properties

% Define the refractive index of each lens 
n1 = getRefractiveIndex(material,lambda);
n2 = getRefractiveIndex(material,lambda);

%% Get remapping function

a1 = 1; % Radius of the input lens 
a2 = 1; % Radius of the output lens

% Gaussian remapping function
sigma = 0.7; % Standard deviation of the Gaussian
[r1,r2] = gaussianRemappingPIAA(a1,a2,sigma,Npts);

%% Make the PIAA sag profiles 

PIAA = makePIAAlenses(r1,r2,n1,n2,L);

%% Ray tracing 

[RAYS,PIAA] = rayTracePIAA(PIAA,Nrays,false);


%% Make Plot

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
    plot([zLens1(end)-lensThickness zLens1(end)],[a1 a1],'k','LineWidth',2);
    plot([zLens1(end)-lensThickness zLens1(end)],[-a1 -a1],'k','LineWidth',2);
    plot([zLens1(end)-lensThickness zLens1(end)-lensThickness],[-a1 a1],'k','LineWidth',2);
    plot([zLens2(end) zLens2(end)+lensThickness],[a2 a2],'k','LineWidth',2);
    plot([zLens2(end) zLens2(end)+lensThickness],[-a2 -a2],'k','LineWidth',2);
    plot([zLens2(end)+lensThickness zLens2(end)+lensThickness],[-a2 a2],'k','LineWidth',2);
    
    % Plot the rays 
    plot(RAYS.z,RAYS.x,'Color',colorOrd(1,:));
    
    xlabel('z/a');
    ylabel('x/a');
    axis equal
    axis([min(RAYS.z(:)) max(RAYS.z(:)) -1.5 1.5])