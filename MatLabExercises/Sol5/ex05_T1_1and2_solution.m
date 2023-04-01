%% Problem Set 5: Solution to 1.1 and 1.2

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

%% 4) Calculate T1S
%s. Eq. (3) on the exercise sheet
    %note that sigma_zz is given in ppm above
    %Hint: Elementwise multiplication in Matlab: array1.*array2

R1S_1 = (deltaIS/4)^2 * (Jiso(omegaI-omegaS, tauC) + 3*Jiso(omegaS, tauC) + 6*Jiso(omegaI+omegaS, tauC)) + 3/4 * (omegaS*sigma_zz1*1e-6).^2 .* Jiso(omegaS, tauC);
R1S_2 = (deltaIS/4)^2 * (Jiso(omegaI-omegaS, tauC) + 3*Jiso(omegaS, tauC) + 6*Jiso(omegaI+omegaS, tauC)) + 3/4 * (omegaS*sigma_zz2*1e-6).^2 .* Jiso(omegaS, tauC);
R1S_3 = (deltaIS/4)^2 * (Jiso(omegaI-omegaS, tauC) + 3*Jiso(omegaS, tauC) + 6*Jiso(omegaI+omegaS, tauC)) + 3/4 * (omegaS*sigma_zz3*1e-6).^2 .* Jiso(omegaS, tauC);

T1S_1 = 1./R1S_1;
T1S_2 = 1./R1S_2;
T1S_3 = 1./R1S_3;

%% 5a) Double Logarithmic Plot of T1S vs omegaS*tauC

figure
    loglog(omegaS*tauC, T1S_1, 'LineWidth', 1.5)
    hold on
    loglog(omegaS*tauC, T1S_2, 'LineWidth', 1.5)
    loglog(omegaS*tauC, T1S_3, 'LineWidth', 1.5)
    xlabel("\omega_S \tau_c")
    ylabel("T_{1, S} [s]")
    set(gca, 'FontSize', 16)
    xlim([tauC(1)*omegaS tauC(end)*omegaS])
    lgd = legend(["0 ppm", "50 ppm", "150 ppm"]);
        lgd.Title.String = "\sigma_{zz}^{(S)}";
        lgd.Location = 'NorthWest';

if(print_to_pdf == 1)
    h = gcf;
    set(h, 'PaperOrientation', 'Landscape')
    set(h, 'PaperUnits', 'normalized')
    set(h, 'PaperPosition', [0 0 1 1])
    print(gcf, '-dpdf', '-painters', './T1S_vs_omegaS_tauC.pdf')
end
    
%% 5a.1) Semi-Log. Plot of R1S vs. tauC for 0 ppm

figure
    semilogx(tauC, R1S_1, 'LineWidth', 1.5)
    xlabel("\tau_c")
    ylabel("R_{1,S} [Hz]")
    xlim([tauC(1) tauC(end)])
    set(gca, 'FontSize', 16)
    
if(print_to_pdf == 1)
    h = gcf;
    set(h, 'PaperOrientation', 'Landscape')
    set(h, 'PaperUnits', 'normalized')
    set(h, 'PaperPosition', [0 0 1 1])
    print(gcf, '-dpdf', '-painters', './R1S_vs_tauC.pdf')
end
    
%% 6) Double Logarithmic Plot of Jiso

J_omegaS = Jiso(omegaS, tauC);
J_omegaSMomegaI = Jiso(omegaS-omegaI, tauC);
J_omegaSPomegaI = Jiso(omegaS+omegaI, tauC);
J_0 = Jiso(0, tauC);

figure
    loglog(omegaS*tauC, J_omegaS, 'LineWidth', 1.5)
    hold on
    loglog(omegaS*tauC, J_omegaSMomegaI, 'LineWidth', 1.5)
    loglog(omegaS*tauC, J_omegaSPomegaI, 'LineWidth', 1.5)
    loglog(omegaS*tauC, J_0, 'LineWidth', 1.5)
    hold off
    xlabel("\omega_S \tau_c")
    ylabel("J (\omega)")
    xlim([tauC(1)*omegaS tauC(end)*omegaS])
    lgd = legend(["\omega_S" ,"\omega_S - \omega_I", "\omega_S + \omega_I", "0"]);
        lgd.Location = 'NorthWest';
    set(gca, 'FontSize', 16)

if(print_to_pdf == 1)
    h = gcf;
    set(h, 'PaperOrientation', 'Landscape')
    set(h, 'PaperUnits', 'normalized')
    set(h, 'PaperPosition', [0 0 1 1])
    print(gcf, '-dpdf', '-painters', './J_vs_tauC.pdf')
end
