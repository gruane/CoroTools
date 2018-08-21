%% Script to compute PIAA lens sag profiles for KPIC
% G. Ruane
% Based on Guyon 2003, Galicher 2005, Vanderbei 2005

clear; close all;
addpath(genpath('PIAA_lib'));

L = 5;% Distance between the PIAA lenses
lambda = 2.2;% Wavlength (microns)
Npts = 10001;% Number of rays 
material = 'CaF2';

label = ['PIAAsag_',material,'_L',num2str(L),'_lam',num2str(lambda),'_Npts',num2str(Npts)];

%% Get material properties

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

% Resample to a single coordinate system for simplicity 
rvals = PIAA.lens1.r;
PIAA.lens2.z = interp1(PIAA.lens2.r,PIAA.lens2.z,rvals,'linear','extrap');

%% Plot the sag profiles 

sag1 = PIAA.lens1.z;
sag2 = PIAA.lens2.z - L;

figure;
plot(rvals,sag1);hold on;
plot(rvals,sag2,'--');hold off;
xlabel('r / a');
ylabel('Sag / a');
legend('Lens 1','Lens 2');

%% Save sag profiles to a file 

M = [rvals',PIAA.lens1.z',(PIAA.lens2.z-L)'];
T = array2table(M,'VariableNames',{'r','z1','z2'});

writetable(T,[label,'.csv']);