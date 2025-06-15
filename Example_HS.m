%% 
% This script demonstrates example usage of the three proposed solvers 
% for finding QSP phase factors. The target function is a polynomial 
% approximation of cos(tau * x), as used in the application of
% Hamiltonian simulation.
%
% The FFPI solver works only in the non-fully coherent setting 
% (i.e., large eta), whereas the other two solvers, HC and INFFT, 
% are applicable in all settings.

addpath("Solvers");
opts.print = 0;
opts.useReal = 0;
opts.targetPre = false;
opts.criteria = 1e-12;
parity = 0;

tau = 100;

%% An fully coherent example
eta = 1e-3;
targ = @(x) (1-eta)*cos(tau.*x);
d = ceil(1.4*tau+log(1e14));
f = chebfun(targ,d);
coef = chebcoeffs(f);
coef = coef(parity+1:2:end)';
    
% Weiss + Half Cholesky algorithm
phi1_HC = HC(coef, eta);

% Weiss + INFFT algorithm
phi1_WeissINFFT = Weiss_INFFT(coef, eta);

%% A non fully coherent example
eta = 0.5;
targ = @(x) (1-eta)*cos(tau.*x);
d = ceil(1.4*tau+log(1e14));
f = chebfun(targ,d);
coef = chebcoeffs(f);
coef = coef(parity+1:2:end)';

% FFPI algorithm
phi2_FFPI = QSP_FFPI(coef',parity,opts);
    
% Weiss + Half Cholesky algorithm
phi2_HC = HC(coef, eta);

% Weiss + INFFT algorithm
phi2_WeissINFFT = Weiss_INFFT(coef, eta);

%% Test whether the results are correct

% For parity = 0, we use the convention that the first reduced phase factor 
% is halved. See details in the paper
% Dong, Yulong, Lin Lin, Hongkang Ni, and Jiasu Wang. "Infinite quantum signal processing." Quantum 8 (2024): 1558.
phi2_HC(1) = phi2_HC(1)/2;
phi2_WeissINFFT(1) = phi2_WeissINFFT(1)/2;

x_grid = -1:0.05:1;
error_FFPI = 0;
error_HC = 0;
error_WeissINFFT = 0;
for x = x_grid
    error_FFPI = max(error_FFPI, abs(ChebyCoef2Func(x, coef, 0, true) - QSPGetPim_sym(phi2_FFPI, x, parity)));
    error_HC = max(error_HC, abs(ChebyCoef2Func(x, coef, 0, true) - QSPGetPim_sym(phi2_HC, x, parity)));
    error_WeissINFFT = max(error_WeissINFFT, abs(ChebyCoef2Func(x, coef, 0, true) - QSPGetPim_sym(phi2_WeissINFFT, x, parity)));
end
disp(['Estimated error of FFPI: ', num2str(error_FFPI)]);
disp(['Estimated error of Weiss + HC: ', num2str(error_HC)]);
disp(['Estimated error of Weiss + INFFT: ', num2str(error_WeissINFFT)]);






