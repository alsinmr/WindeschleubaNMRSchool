%% Problem Set 5: Relaxation Rate Constants: T1
% Problems 1.1 and 1.2

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

%% 4) Calculate T1S
%s. Eq. (3) on the exercise sheet
    %note that sigma_zz is given in ppm above
    %Hint: Elementwise multiplication in Matlab: array1.*array2

R1S = ;

T1S = 1./R1S;
%% 5a) Double Logarithmic Plot of T1S vs omegaS*tauC

figure
    loglog(omegaS*tauC, T1S)
    xlabel("\omega_S \tau_c")
    ylabel("T_{1, S} [s]")
    xlim([tauC(1)*omegaS tauC(end)*omegaS])
    
%% 5a.1) Semi-Log. Plot of R1S vs. tauC

figure
    semilogx(tauC, R1S)
    xlabel("\tau_c")
    ylabel("R_{1,S} [Hz]")
    xlim([tauC(1) tauC(end)])
    
%% 6) Double Logarithmic Plot of Jiso

J_omegaS = ;                      %spectral density at omegaS
J_omegaSMomegaI = ;        %spectral density at omegaS-omegaI
J_omegaSPomegaI = ;        %spectral density at omegaS+omegaI
J_0 = ;                                %spectral density at 0

figure
    loglog(omegaS*tauC, J_omegaS)
    hold on
    loglog(omegaS*tauC, J_omegaSMomegaI)
    loglog(omegaS*tauC, J_omegaSPomegaI)
    loglog(omegaS*tauC, J_0)
    hold off
    xlabel("\omega_S \tau_c")
    ylabel("J (\omega)")
    legend(["\omega_S" ,"\omega_S - \omega_I", "\omega_S + \omega_I", "0"])


