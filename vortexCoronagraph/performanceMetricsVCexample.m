%% Computes performance metrics for a vortex coronagraph
% author: G. Ruane 

clear; close all;
addpath('utils');

N = 2^12;% Size of the computational grid (aperture will be padded to NxN)

lambdaOverD = 4; % lambda/D in focal plane (units of samples)
apRad = N/2/lambdaOverD; % Aperture radius in units of samples 

R0 = 0.1; % Central obscuration radius (units of outer radius)
LSin = 0.3; % Inner radius of the Lyot stop 
LSout = 0.95; % Outer radius of the Lyot stop

Qin = 2;  % Inner radius of image plane "dark hole" (units of lambdaOverD)
Qout= 10; % Outer radius of image plane "dark hole" (units of lambdaOverD)

radiusOfPSFROI = 1*lambdaOverD; % Choose a radius for the PSF "region of interest" 
c = 50; % hypergaussian ROI power 

charge = 4; % Charge of the vortex focal plane mask

useGPU = false; % Use the GPU? Keep this false unless you know what you're doing 

%% Initialize variables

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

iPSF_unnormalized = abs(myfft2(EP)).^2; % PSF without normalization 
normI = max(iPSF_unnormalized(:));% Irradiance normalization factor
totEnergy = sum(iPSF_unnormalized(:)); % Total energy in telescope PSF

iPSFtel = iPSF_unnormalized/normI; % Normalized irradiance in the image plane 

Q = and(RHO > Qin*lambdaOverD, RHO < Qout*lambdaOverD); % "Dark hole" region mask 

PSFROI = exp(-(RHO/radiusOfPSFROI).^c); % "Region of interest" or "photometric aperture" 

FPM = exp(1i*charge*THETA); % Focal plane mask

% Read in optimal parameters for propagating through the vortex
if(abs(charge) > 0)
    [ inVal, outVal ] = loadOptVals( charge, N );
else
    inVal = 0.1;
    outVal= 1.22;
end

% Show pupil field
figure;
    imagesc(xvals/apRad,yvals/apRad,abs(EP)+abs(LS));
    axis image;
    axis([-1 1 -1 1]);
    title('Pupil and Lyot stop');
    hx = xlabel('{\itx} / {\itR}');
    hy = ylabel('{\ity} / {\itR}');
    set(gca,'XTick',-1:0.5:1,'YTick',-1:0.5:1);
    set(gca,'TickDir','out','YDir','normal');
    drawnow;
    

%% Compute eta_s

% Step 1: Compute the on-axis PSF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes the Lyot plane field 
LP = vortexCoronagraph_Pup2Pup( EP, FPM, apRad, lambdaOverD, RHO, N, 'dft', 'forward', inVal, outVal, useGPU );

FP = myfft2(LP.*LS); % Computes the focal plane field 

iPSF = abs(FP).^2/normI; % Normalized irradiance in the image plane 

% Show PSF
figure;
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(iPSF));
    colorbar; 
    if(max(log10(iPSF(:)))>-10)
        caxis([-10 max(log10(iPSF(:)))]);
    end
    axis image;
    axis([-20 20 -20 20]);
    title('Stellar PSF, normalized irradiance (log scale)');
    hx = xlabel('Angular coordinate (\lambda/D)');
    hy = ylabel('Angular coordinate (\lambda/D)');
    set(gca,'XTick',-20:5:20,'YTick',-20:5:20);
    set(gca,'TickDir','out','YDir','normal');
    drawnow;

% Step 2: Convolve with the PSFROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta_s = myifft2(myfft2(iPSF*normI/totEnergy).*myfft2(PSFROI));

% Show convolved PSF 
figure;
    imagesc(xvals/lambdaOverD,yvals/lambdaOverD,log10(eta_s));
    colorbar;
	if(max(log10(iPSF(:)))>-10)
        caxis([-10 max(log10(iPSF(:)))]);
    end
    axis image;
    axis([-20 20 -20 20]);
    title('Fraction of starlight detected in ROI ( \eta_s )');
    hx = xlabel('Angular coordinate (\lambda/D)');
    hy = ylabel('Angular coordinate (\lambda/D)');
    set(gca,'XTick',-20:5:20,'YTick',-20:5:20);
    set(gca,'TickDir','out','YDir','normal');
    drawnow;

mean_eta_s = mean(eta_s(Q)) % Computes mean eta_s within the "dark hole" 

%% Compute eta_p (in 1D)

offsets = (0:20)*lambdaOverD; % Offsets in the units of image plane samples

eta_p = []; % Allocate array to hold computed eta_p values 
for offset = offsets
    
    tilt = exp(1i*2*pi*offset/N*X); % Tilted wavefront at pupil
    PSFROIoffaxis = exp(-(sqrt((X-offset).^2+Y.^2)/radiusOfPSFROI).^c); % Shifted PSF "region of interest" 

    % The inner most region is computed with oversampled dfts instead of ffts 
	if(offset/lambdaOverD < 1)
        method = 'dft';
    else
        method = 'fft';
    end
    if(useGPU)
        FP = gather(FP);
    end
    
    % Compute PSF
    LP = vortexCoronagraph_Pup2Pup( EP.*tilt, FPM, apRad, lambdaOverD, RHO, N, method, 'forward', inVal, outVal, useGPU );
    FP = myfft2(LP.*LS); % Computes the focal plane field 
    iPSFoffaxis = abs(FP).^2/normI; % Normalized irradiance in the image plane 
    
    % Plot off-axis PSFs while calculating the eta_p
    figure(901);
        subplot(1,2,1);
        imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFoffaxis);
        caxis([0 1]);
        axis image;
        axis([-20 20 -20 20]);
        title('"Planet" light');
        hx = xlabel('Angular coordinate (\lambda/D)');
        hy = ylabel('Angular coordinate (\lambda/D)');
        set(gca,'XTick',-20:5:20,'YTick',-20:5:20);
        set(gca,'TickDir','out','YDir','normal');
        
        subplot(1,2,2);
        imagesc(xvals/lambdaOverD,yvals/lambdaOverD,iPSFoffaxis.*PSFROIoffaxis);
        caxis([0 1]);
        axis image;
        axis([-20 20 -20 20]);
        title('"Planet" light inside the ROI');
        hx = xlabel('Angular coordinate (\lambda/D)');
        hy = ylabel('Angular coordinate (\lambda/D)');
        set(gca,'XTick',-20:5:20,'YTick',-20:5:20);
        set(gca,'TickDir','out','YDir','normal');
        drawnow;
    
    % Calculate eta_p and append it to the eta_p array 
    eta_p = [eta_p,sum(sum(PSFROIoffaxis.*abs(FP).^2))/totEnergy];
        
end

% Maximum eta_p with a circular pupil 
eta_p_ideal = (LSout)^2*(1-besselj(0,radiusOfPSFROI/lambdaOverD*pi*LSout)^2-besselj(1,radiusOfPSFROI/lambdaOverD*pi*LSout)^2);

% Plot eta_p
figure;
    plot(offsets/lambdaOverD,eta_p,'-o'); hold on;
    plot([0 offsets(end)/lambdaOverD],[eta_p_ideal eta_p_ideal],'--'); hold off; 
    xlabel('Angular separationQ (\lambda/D)');
    ylabel('\eta_p (along horizontal axis)');
    title(['Absolute throughput within ',num2str(radiusOfPSFROI/lambdaOverD),' \lambda/D radius of the planet position']);
    ylim([0 1]);
    
%% Raw contrast and relative exposure time

% Compute eta_p and eta_s for the telescope 
eta_p_tel = sum(sum(PSFROI.*iPSF_unnormalized))/totEnergy;
eta_s_tel = myifft2(myfft2(iPSFtel*normI/totEnergy).*myfft2(PSFROI));

% Take the azimuthal average of eta_s in each case (can also use std)
eta_s_1D     = mean(polarTransform(eta_s    ,1,1),2)';
eta_s_tel_1D = mean(polarTransform(eta_s_tel,1,1),2)';

% Interpolate eta_p
eta_p_interp = interp1(offsets,eta_p,0:N/2-1,'spline','extrap');

% Compute the relative exposure time for detection in the
% photon-noise-limited regime 
relExpTime1D = (eta_s_1D./eta_p_interp.^2)./(eta_s_tel_1D./eta_p_tel.^2);

figure;
    semilogy((0:N/2-1)/lambdaOverD,eta_s_tel_1D,'-o');hold on;
    semilogy((0:N/2-1)/lambdaOverD,eta_s_1D,'-o');hold off;
    xlabel('Angular separation (\lambda/D)');
    ylabel('\eta_s');
    legend('w/o coro','w/ coro');
    xlim([0 20]);

figure;
    semilogy((0:N/2-1)/lambdaOverD,eta_s_tel_1D/eta_p_tel);hold on;
    semilogy((0:N/2-1)/lambdaOverD,eta_s_1D./eta_p_interp);hold off;
    xlabel('Angular separation (\lambda/D)');
    ylabel('Raw contrast: \eta_s / \eta_p');
    title('Raw contrast');
    legend('w/o coro','w/ coro');
    xlim([0 20]);

figure;
    semilogy((0:N/2-1)/lambdaOverD,relExpTime1D);
    xlabel('Angular separation (\lambda/D)');
    ylabel('\Deltat_{coro} / \Deltat_{tel}');
    title('Relative exposure time');
    xlim([1 20]);

