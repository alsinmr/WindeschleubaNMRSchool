function T1x = T1x(delta, omegaI, omegaS, tau)
% calculates T1x
T1x = ((delta/4)^2*(-1*J(omegaI-omegaS,tau) + ... 
6*J(omegaI+omegaS,tau))).^-1;
