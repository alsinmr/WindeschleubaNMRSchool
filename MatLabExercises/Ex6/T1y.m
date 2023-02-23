function T1y = T1y(delta, sigmazzS, omegaS, tau)
% calculates T1 time 
T1y = (3/4*delta*sigmazzS*omegaS*J(omegaS,tau)).^-1;