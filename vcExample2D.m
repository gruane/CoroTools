%% Computes the amount of leaked starlight through a vortex coronagraph 
% author: G. Ruane 

clear;

addpath('utils');

N = 2^12;% Size of the computational grid (aperture will be padded to NxN)

lambdaOverD = 4; % lambda/D in focal plane (units of samples)
apRad = N/2/lambdaOverD; % Aperture radius in units of samples 

R0 = 0.1; % Central obscuration radius (units of outer radius)
LSin = 0.2; % Inner radius of the Lyot stop 
LSout = 0.95; % Outer radius of the Lyot stop

Qin = 2.5; % Inner radius of image plane region-of-interest (units of lambdaOverD)
Qout= 3.5; % Outer radius of image plane region-of-interest (units of lambdaOverD)

charge = 4; % Charge of the focal plane mask

useGPU = false; % Use the GPU? Keep this false unless you know what you're doing 

%%

% Defines the coordinate systems
[X,Y] = meshgrid(-N/2:N/2-1); % Grids with Cartesian (x,y) coordinates 
[THETA,RHO] = cart2pol(X,Y);  % Grids with polar (rho,theta) coordinates 
xvals = X(1,:);yvals = Y(:,1);

EP = exp(-(RHO/(apRad)).^1000); % Entrance pupil function
LS = exp(-(RHO/(LSout*apRad)).^1000); % Lyot stop function 

% Add central obscuration, if neccessary 
if(R0 > 0)
    EP = EP - exp(-(RHO/(R0*apRad)).^1000); % Entrance pupil function w. central obscuration
end
if(LSin > 0)
    LS = LS - exp(-(RHO/(LSin*apRad)).^1000); % Lyot stop function w. central obscuration
end

normI = max(max(abs(myfft2(EP)).^2));% Irradiance normalization factor

Q = and(RHO > Qin*lambdaOverD, RHO < Qout*lambdaOverD); % ROI mask 
    
FPM = exp(1i*charge*THETA); % Focal plane mask

tic;
if(abs(charge) > 0)
    [ inVal, outVal ] = loadOptVals( charge, N );
else
    inVal = 0.1;
    outVal= 1.22;
end
% Computes the Lyot plane mode 
LP = vortexCoronagraph_Pup2Pup( EP, FPM, apRad, lambdaOverD, RHO, N, 'dft', 'forward', inVal, outVal, useGPU );
toc; 

FP = myfft2(LP.*LS); % Computes the focal plane field 

iPSF = abs(FP).^2/normI; % Normalized irradiance in the image plane 

irradianceROI = mean(iPSF(Q)) % Computes the normalized irradiance at the inner working angle 

%% Plot the result  

figure(1); 

subplot(1,2,1);
imagesc(xvals/apRad,yvals/apRad,abs(EP)+abs(LS));
colorbar; 
axis image;
axis([-1 1 -1 1]);
title('Pupil and Lyot stop');
hx = xlabel('{\itx} / {\itR}');
hy = ylabel('{\ity} / {\itR}');
set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
set(gca,'TickDir','out');set(gca,'YDir','normal');


subplot(1,2,2);
imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(iPSF));
colorbar; caxis([-10 max(log10(iPSF(:)))]);
axis image;
axis([-10 10 -10 10]);
title('Normalized irradiance (log scale)');
hx = xlabel('Angular coordinate (\lambda/D)');
hy = ylabel('Angular coordinate (\lambda/D)');
set(gca,'XTick',-10:5:10,'YTick',-10:5:10);
set(gca,'TickDir','out');set(gca,'YDir','normal');
