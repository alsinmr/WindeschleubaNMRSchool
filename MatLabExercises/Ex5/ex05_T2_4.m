%% Problem Set 5: Backcalculate Correlation Times
%Problem 2.4

%% 0) Clear Workspace and Close Figures
clear, close all

%% 1) Define Constants

mu0 = 4*pi*1e-7;                    %V s / A / m
hbar = 1.05457266e-34;              %J s / rad
h = 2*pi*hbar;                      %J s
gammaH = 26.7522128e7;              %rad / s / T
gammaN = -2.7116e7;                 %rad / s / T
rNH = 1.1*1e-10;                    %m 

%% 2) Functions

%Spectral density function (for isotropic tumbling)
    %Call as: Jiso(omega, tau)
    %one of the two inputs can be an array
Jiso = @(omega, tau) (2/5)*tau./(1+(omega.*tau).^2);

%% 3) Parameters

B0 = ;                          %T
sigma_zz = ;                       %ppm
tauC = ;           %s

omegaI = ;                 %rad/s, 1H Larmor frequency
omegaS = ;                 %rad/s, 15N Larmor frequency

deltaIS = -2*mu0*gammaH*gammaN*hbar/(4*pi)/rNH^3;

%% 4) Calculate T1S and T2S
%s. Eqs. (3) and (4) on the exercise sheet
    %note that sigma_zz is given in ppm above
    %Hint: Elementwise multiplication in Matlab: array1.*array2

R1S = ;

T1S = 1./R1S;

R2S = ;

T2S = 1./R2S;

%% 5) Backcalculate tauC
%s. Eq. (5) on the exercise sheet
    %Hint: Elementwise division in Matlab: array1./array2
    
tauC_BC = ;

%% 6) Double Logarithmic Plot of backcalculated tauC against input tauC

figure
    loglog(tauC, real(tauC_BC))
    hold on
    loglog(tauC, tauC)
    xlabel("\tau_c [s]")
    ylabel("\tau_c [s]")
    lgd = legend(["\tau_c (input) vs. \tau_c (fit)", "\tau_c (input) vs. \tau_c (input)"]);
    xlim([tauC(1) tauC(end)])