function [ inVal, outVal ] = loadOptVals( charge, N )
%loadOptVals Returns the optimal window values for a given array size. 
%   Helper function for vortexCoronagraph_Pup2Pup
%   Inputs: 'charge' - Charge of the vortex
%           'N' - Size of the padded array 

    Ns = 2.^(10:13);
    load('optVals.mat');
    row = charge/2;
    dim = round(interp1(Ns,1:numel(Ns),N));
    inVal = M(row,1,dim);
    outVal = M(row,2,dim);
    
end

