%% Problem Set 5: Solution to 2.4

%% 0) Clear Workspace and Close Figures
clear, close all

%% 1) Define Constants

mu0 = 4*pi*1e-7;                    %V s / A / m
hbar = 1.05457266e-34;              %J s / rad
h = 2*pi*hbar;                      %J s
gammaH = 26.7522128e7;              %rad / s / T
gammaN = -2.7116e7;                 %rad / s / T
rNH = 1.1*1e-10;                    %m 

print_to_pdf = 0;                   %set to 1 to print plots to pdf

%% 2) Functions

%Spectral density function (for isotropic tumbling)
    %Call as: Jiso(omega, tau)
    %one of the two inputs can be an array
Jiso = @(omega, tau) (2/5)*tau./(1+(omega.*tau).^2);

%% 3) Parameters

B0 = 14.1;                          %T
sigma_zz = 0;                       %ppm
tauC = 10.^(-12:0.01:-5);           %s

omegaI = -B0*gammaH;                 %rad/s, 1H Larmor frequency
omegaS = -B0*gammaN;                 %rad/s, 15N Larmor frequency

deltaIS = -2*mu0*gammaH*gammaN*hbar/(4*pi)/rNH^3;

%% 4) Calculate T1S and T2S
%s. Eqs. (3) and (4) on the exercise sheet
    %note that sigma_zz is given in ppm above
    %Hint: Elementwise multiplication in Matlab: array1.*array2

R1S = (deltaIS/4)^2 * (Jiso(omegaI-omegaS, tauC) + 3*Jiso(omegaS, tauC) + 6*Jiso(omegaI+omegaS, tauC)) + 3/4 * (omegaS*sigma_zz*1e-6).^2 .* Jiso(omegaS, tauC);

T1S = 1./R1S;

R2S = deltaIS^2/32 * (4*Jiso(0,tauC) + Jiso(omegaI-omegaS,tauC) + 3*Jiso(omegaS,tauC) + 6*Jiso(omegaI,tauC) + 6*Jiso(omegaI+omegaS,tauC))...
       + 1/2 *(omegaS*sigma_zz*1e-6).^2 *Jiso(0,tauC) + 3/8 *(omegaS*sigma_zz*1e-6).^2 * Jiso(omegaS,tauC);

T2S = 1./R2S;

%% 5) Backcalculate tauC
%s. Eq. (5) on the exercise sheet
    %Hint: Elementwise division in Matlab: array1./array2
    
tauC_BC = 1/2 * sqrt((1/omegaS^2) *(6*T1S./T2S - 7));

%% 6) Double Logarithmic Plot of backcalculated tauC against input tauC

figure
    loglog(tauC, real(tauC_BC), 'LineWidth', 1.5)
    hold on
    loglog(tauC, tauC, 'LineWidth', 1.5)
    xlabel("\tau_c [s]")
    ylabel("\tau_c [s]")
    set(gca, 'FontSize', 16)
    lgd = legend(["\tau_c (input) vs. \tau_c (fit)", "\tau_c (input) vs. \tau_c (input)"]);
    xlim([tauC(1) tauC(end)])
    
if(print_to_pdf == 1)
    h = gcf;
    set(h, 'PaperOrientation', 'Landscape')
    set(h, 'PaperUnits', 'normalized')
    set(h, 'PaperPosition', [0 0 1 1])
    print(gcf, '-dpdf', '-painters', './tauC_backcalculated.pdf')
end    