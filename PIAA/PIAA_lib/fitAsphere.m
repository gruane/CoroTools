function [sag1,sag2,PIAA] = fitAsphere(PIAA,polyOrders,rvals)
%[sag1,sag2,PIAA] = fitAsphere(PIAA,polyOrders,rvals)
%   Fits an asphere profile of order 'polyOrders'. Returns the sag profiles
%   and puts the best fit parameters into the PIAA struct
%
%   Inputs:
%       PIAA - Structure containing the PIAA parameters from makePIAAlenses.m
%       polyOrders - Number of polynomial terms
%       rvals - Radius values to evaluate the fit 
%
%   Outputs:
%       sag1, sag1 - Sag of lens 1 and 2, respectively

    X0 = [-10,0.01,0, zeros(1,polyOrders)];% Initial guess 
    
    options = optimoptions('lsqcurvefit','Display','off');
    % Least squares fit to asphereEqn
    PIAA.lens1.asphFitParams = lsqcurvefit(@asphereEqn,X0,PIAA.lens1.r,PIAA.lens1.z,[],[],options);
    PIAA.lens2.asphFitParams = lsqcurvefit(@asphereEqn,X0,PIAA.lens2.r,PIAA.lens2.z - PIAA.L,[],[],options);

    % Evaluate sag fits at rvals
    sag1 = asphereEqn(PIAA.lens1.asphFitParams,rvals);
    sag2 = asphereEqn(PIAA.lens2.asphFitParams,rvals);
    
end

