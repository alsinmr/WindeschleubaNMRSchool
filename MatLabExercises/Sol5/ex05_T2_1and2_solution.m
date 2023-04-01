%% Problem Set 5: Solution to 2.1 and 2.2

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
sigma_zz1 = 0;                      %ppm
sigma_zz2 = 50;                     %ppm
sigma_zz3 = 150;                    %ppm
tauC = 10.^(-12:0.01:-5);           %s

omegaI = -B0*gammaH;                 %rad/s, 1H Larmor frequency
omegaS = -B0*gammaN;                 %rad/s, 15N Larmor frequency

deltaIS = -2*mu0*gammaH*gammaN*hbar/(4*pi)/rNH^3;

%% 4) Calculate T2S
%s. Eq. (4) on the exercise sheet
    %note that sigma_zz is given in ppm above
    %Hint: Elementwise multiplication in Matlab: array1.*array2

R2S_1 = deltaIS^2/32 * (4*Jiso(0,tauC) + Jiso(omegaI-omegaS,tauC) + 3*Jiso(omegaS,tauC) + 6*Jiso(omegaI,tauC) + 6*Jiso(omegaI+omegaS,tauC))...
       + 1/2 *(omegaS*sigma_zz1*1e-6).^2 *Jiso(0,tauC) + 3/8 *(omegaS*sigma_zz1*1e-6).^2 * Jiso(omegaS,tauC);
R2S_2 = deltaIS^2/32 * (4*Jiso(0,tauC) + Jiso(omegaI-omegaS,tauC) + 3*Jiso(omegaS,tauC) + 6*Jiso(omegaI,tauC) + 6*Jiso(omegaI+omegaS,tauC))...
       + 1/2 *(omegaS*sigma_zz2*1e-6).^2 *Jiso(0,tauC) + 3/8 *(omegaS*sigma_zz2*1e-6).^2 * Jiso(omegaS,tauC);
R2S_3 = deltaIS^2/32 * (4*Jiso(0,tauC) + Jiso(omegaI-omegaS,tauC) + 3*Jiso(omegaS,tauC) + 6*Jiso(omegaI,tauC) + 6*Jiso(omegaI+omegaS,tauC))...
       + 1/2 *(omegaS*sigma_zz3*1e-6).^2 *Jiso(0,tauC) + 3/8 *(omegaS*sigma_zz3*1e-6).^2 * Jiso(omegaS,tauC);

   
T2S_1 = 1./R2S_1;
T2S_2 = 1./R2S_2; 
T2S_3 = 1./R2S_3;

%% 5) Double Logarithmic Plot of T2S

figure
    loglog(omegaS*tauC, T2S_1, 'LineWidth', 1.5)
    hold on
    loglog(omegaS*tauC, T2S_2, 'LineWidth', 1.5)
    loglog(omegaS*tauC, T2S_3, 'LineWidth', 1.5)
    xlabel("\omega_S \tau_c")
    ylabel("T_{2, S} [s]")
    xlim([omegaS*tauC(1) omegaS*tauC(end)])
    set(gca, 'FontSize', 16)
    lgd = legend(["0 ppm", "50 ppm", "150 ppm"]);
        lgd.Title.String = "\sigma_{zz}^{(S)}";
        lgd.Location = 'Northeast';
        
if(print_to_pdf == 1)
    h = gcf;
    set(h, 'PaperOrientation', 'Landscape')
    set(h, 'PaperUnits', 'normalized')
    set(h, 'PaperPosition', [0 0 1 1])
    print(gcf, '-dpdf', '-painters', './T2S_vs_omegaS_tauC.pdf')
end
