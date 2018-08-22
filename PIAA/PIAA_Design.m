%% Script to compute PIAA lens sag profiles for KPIC
% G. Ruane
% Based on Guyon 2003, Galicher 2005, Vanderbei 2005

clear; close all;
addpath(genpath('PIAA_lib'));
addpath(genpath('export_scripts'));

L = 10;% Distance between the PIAA lenses
lambda = 2.2;% Wavlength (microns)
Npts = 10001;% Number of rays for design 
Npts_resamp = 20001;% Number of rays for final resampling  
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

%% Fit asphere equation 

polyOrders = 16; % Number of asphere polynomial terms (a_4*R^4 + a_6*R^6 + ...)
rvals = linspace(0, 1, Npts_resamp);% r values to evaluate fit 
[sagFitlens1,sagFitlens2,PIAA] = fitAsphere(PIAA,polyOrders,rvals);

spherRef1 = asphereEqn(PIAA.lens1.asphFitParams(1:3),rvals);
spherRef2 = asphereEqn(PIAA.lens2.asphFitParams(1:3),rvals);

%% Plot the sag profiles and fits 

% plotPIAAdesign;

figure;
plot(PIAA.lens1.r,PIAA.lens1.z,'LineWidth',2); hold on;
plot(PIAA.lens2.r,PIAA.lens2.z - L,'--','LineWidth',2);
plot(rvals,sagFitlens1);
plot(rvals,sagFitlens2);
plot(rvals,spherRef1,':');
plot(rvals,spherRef2,':');
hold off;
xlabel('r / a');
ylabel('Sag / a');
legend('Lens 1','Lens 2');
axis([0 1 min([sagFitlens1 sagFitlens2]) max([sagFitlens1 sagFitlens2])]);

%% Save sag profiles to a file 

% Resample to a single coordinate system for simplicity 
sag1 = interp1(PIAA.lens1.r,PIAA.lens1.z,rvals,'linear','extrap');
sag2 = interp1(PIAA.lens2.r,PIAA.lens2.z - L,rvals,'linear','extrap');

M = [rvals',sag1',sag2',sagFitlens1',sagFitlens2'];
T = array2table(M,'VariableNames',{'r','sag1','sag2','asphFit1','asphFit2'});

writetable(T,[label,'.csv']);