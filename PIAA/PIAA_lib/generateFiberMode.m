function mode = generateFiberMode( diam_samples, coords )
%mode = generateFiberMode( diam_samples, coords )
%   Generates the Gaussian fiber mode function

mode = sqrt(2/(pi*(diam_samples/2)^2))* ...% Normalization factor 
            exp(-(coords.RHO/(diam_samples/2)).^2);% Gaussian 

end

