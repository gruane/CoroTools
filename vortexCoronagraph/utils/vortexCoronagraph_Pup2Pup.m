function OUT = vortexCoronagraph_Pup2Pup( IN, FPM, apRad, lambdaOverD, RHO, N, algo, operation, inVal, outVal, useGPU )
%vortexCoronagraph_Pup2Pup Propagates from the input pupil to output pupil
%   Inputs: 
%       'IN' - 2D array of complex field values in the input pupil.
%       'FPM' - 2D complex representation of the vortex focal plane mask
%           exp(1i*l*theta).
%       'apRad' - Number of samples across the aperture radius.
%       'lambdaOverD' - Number of samples per lambda/D after an FFT 
%       'RHO' - 2D array of radial coordinates (units of samples) 
%       'N' - Size of the padded array IN
%       'algo' - string - 'fft' or 'dft' 
%       'operation' - 'fwd' or 'adj' for the forward or adjoint operation 
%       'inVal' and 'outVal' define regions of the image plane with fine
%           and coarse sampling. loadOptVals returns "optimal" ones for select
%           array sizes. 
%       'useGPU' - Runs demanding calculations on the GPU (as Matlab gpuArrays)


    showPlots2debug = false; 

    if(useGPU)
        IN = gpuArray(IN);
    end
    if( strcmp(algo,'fft') && ~strcmp(operation,'adj') )
        
        EP = IN;
        OUT = myifft2(myfft2(EP).*FPM);
        
	elseif( strcmp(algo,'fft') && strcmp(operation,'adj') )
        
        LP = IN;
        OUT = myifft2(myfft2(LP).*conj(FPM));
        
    elseif( strcmp(algo,'dft') )

        cut_rad1 = inVal*lambdaOverD;
        cut_rad2 = outVal*lambdaOverD;
        
        windowKnee = 1-cut_rad1/cut_rad2;

        D = 2*apRad;
        
        NA = 2^(nextpow2(D))+2;
        crop = N/2-NA/2+1:N/2+NA/2;
        
        
        % DFT vectors 
        x = ((0:NA-1)-NA/2)/D;
        u1 = ((0:N-1)-N/2)/lambdaOverD;
        u2 = ((0:N-1)-N/2)*2*cut_rad2/N;
        
        windowMASK1 = generateTukeyWindow( 2*cut_rad2*lambdaOverD, RHO, windowKnee ) ;
        windowMASK2 = generateTukeyWindow( N, RHO, windowKnee ) ;
        
        if(useGPU)
            x = gpuArray(x);
            u1 = gpuArray(u1);
            u2 = gpuArray(u2);
            windowMASK1 = gpuArray(windowMASK1);
            windowMASK2 = gpuArray(windowMASK2);
        end
        if(~strcmp(operation,'adj'))

            EP = IN;
            EP = EP(crop,crop);
            if showPlots2debug; figure;imagesc(abs(EP));axis image;colorbar; title('Cropped pupil'); end;
            
            %%%%%%% Large scale DFT

            FP1 = (N/lambdaOverD)/(D*N)*exp(-1i*2*pi*u1'*x)*EP*exp(-1i*2*pi*x'*u1); 
            if showPlots2debug; figure;imagesc(log10(abs(FP1).^2));axis image;colorbar; title('Large scale DFT'); end;
            LP1 = (N/lambdaOverD)/(D*N)*exp(1i*2*pi*x'*u1)*(FP1.*FPM.*(1-windowMASK1))*exp(1i*2*pi*u1'*x);

            %%%%%%% Fine sampled DFT

            FP2 = 2*cut_rad2/(D*N)*exp(-1i*2*pi*u2'*x)*EP*exp(-1i*2*pi*x'*u2); 
            if showPlots2debug; figure;imagesc(log10(abs(FP2).^2));axis image;colorbar; title('Fine sampled DFT'); end;
            LP2 = 2*cut_rad2/(D*N)*exp(1i*2*pi*x'*u2)*(FP2.*FPM.*windowMASK2)*exp(1i*2*pi*u2'*x);

            OUT = padarray_centered(LP1+LP2,NA,NA,N);
            %disp('Propagating through vortex with forward DFT.');
            if showPlots2debug; figure;imagesc(abs(LP1+LP2));axis image;colorbar; title('Lyot plane'); end;
            if showPlots2debug; figure;imagesc(abs(LP1+LP2-EP));axis image;colorbar; title('Lyot plane - Entrance Pupil'); end;
        elseif(strcmp(operation,'adj') )
        
            LP = IN(crop,crop);

            %%%%%%% Large scale DFT

            FP1 = (N/lambdaOverD)/(D*N)*exp(-1i*2*pi*u1'*x)*LP*exp(-1i*2*pi*x'*u1); 
            EP1 = (N/lambdaOverD)/(D*N)*exp(1i*2*pi*x'*u1)*(FP1.*conj(FPM).*(1-windowMASK1))*exp(1i*2*pi*u1'*x);


            %%%%%%% Fine sampled DFT

            FP2 = 2*cut_rad2/(D*N)*exp(-1i*2*pi*u2'*x)*LP*exp(-1i*2*pi*x'*u2); 
            EP2 = 2*cut_rad2/(D*N)*exp(1i*2*pi*x'*u2)*(FP2.*conj(FPM).*windowMASK2)*exp(1i*2*pi*u2'*x);

            OUT = padarray_centered(EP1+EP2,NA,NA,N);
        end
    else
        error('Error. \nChoose algo = fft or dft. \nChoose forward or adj operation.')
    end
    if(useGPU)
        OUT = gather(OUT);
    end
end

