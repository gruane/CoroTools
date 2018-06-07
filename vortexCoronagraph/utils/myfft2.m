function out = myfft2( in )
%myfft2 Performs 2D FFT with fftshifts before and after calling fft2.

    out = fftshift(fft2(ifftshift(in)));

end

