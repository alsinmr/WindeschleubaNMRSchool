function Jiso = J(omega, tau)
% calculates Jiso
Jiso = 2/5 * tau ./(1+(omega .* tau).^2);