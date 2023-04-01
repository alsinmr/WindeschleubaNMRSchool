%% Problem Set 5: Solution to 2.3

%% 0) Clear Workspace and Close Figures
clear, close all

%% 1) Define Constants

mu0 = 4*pi*1e-7;                    %V s / A / m
hbar = 1.05457266e-34;              %J s / rad
h = 2*pi*hbar;                      %J s
gammaH = 26.7522128e7;              %rad / s / T
gammaN = -2.7116e7;                 %rad / s / T
rNH = 1.1*1e-10;                    %m 

print_to_pdf = 0;                   %set to one to print plots to pdf

%% 2) Functions

%Spectral density function (for isotropic tumbling)
    %Call as: Jiso(omega, tau)
    %one of the two inputs can be an array
Jiso = @(omega, tau) (2/5)*tau./(1+(omega.*tau).^2);

%% 3) Parameters

B0 = 1:0.1:30;                       %T
sigma_zz = 0;                        %ppm
tauC1 = 10^-12;                      %s
tauC2 = 10^-10;                      %s
tauC3 = 10^-8;                       %s

omegaI = -B0*gammaH;                 %rad/s, 1H Larmor frequency
omegaS = -B0*gammaN;                 %rad/s, 15N Larmor frequency

deltaIS = -2*mu0*gammaH*gammaN*hbar/(4*pi)/rNH^3;

%% 4) Calculate T2S
%s. Eq. (4) on the exercise sheet
    %note that sigma_zz is given in ppm above
    %Hint: Elementwise multiplication in Matlab: array1.*array2
    
R2S_tauC1 = deltaIS^2/32 * (4*Jiso(0,tauC1) + Jiso(omegaI-omegaS,tauC1) + 3*Jiso(omegaS,tauC1) + 6*Jiso(omegaI,tauC1) + 6*Jiso(omegaI+omegaS,tauC1))...
       + 1/2 *(omegaS*sigma_zz*1e-6).^2 *Jiso(0,tauC1) + 3/8 *(omegaS*sigma_zz*1e-6).^2 .* Jiso(omegaS,tauC1);
R2S_tauC2 = deltaIS^2/32 * (4*Jiso(0,tauC2) + Jiso(omegaI-omegaS,tauC2) + 3*Jiso(omegaS,tauC2) + 6*Jiso(omegaI,tauC2) + 6*Jiso(omegaI+omegaS,tauC2))...
       + 1/2 *(omegaS*sigma_zz*1e-6).^2 *Jiso(0,tauC2) + 3/8 *(omegaS*sigma_zz*1e-6).^2 .* Jiso(omegaS,tauC2);
R2S_tauC3 = deltaIS^2/32 * (4*Jiso(0,tauC3) + Jiso(omegaI-omegaS,tauC3) + 3*Jiso(omegaS,tauC3) + 6*Jiso(omegaI,tauC3) + 6*Jiso(omegaI+omegaS,tauC3))...
       + 1/2 *(omegaS*sigma_zz*1e-6).^2 *Jiso(0,tauC3) + 3/8 *(omegaS*sigma_zz*1e-6).^2 .* Jiso(omegaS,tauC3);

T2S_tauC1 = 1./R2S_tauC1;
T2S_tauC2 = 1./R2S_tauC2;
T2S_tauC3 = 1./R2S_tauC3;
   
%% 5 Semi-Logarithmic Plot of T2S vs B0

figure
    semilogy(B0, T2S_tauC1, 'LineWidth', 1.5)
    hold on
    semilogy(B0, T2S_tauC2, 'LineWidth', 1.5)
    semilogy(B0, T2S_tauC3, 'LineWidth', 1.5)
    xlabel("B_0 [T]")
    ylabel("T_{2, S} [s]")
    lgd = legend([num2str(tauC1); num2str(tauC2); num2str(tauC3)]);
        lgd.Title.String = "\tau_c";
        lgd.Location = 'Northeast';
    xlim([B0(1) B0(end)])
    ylim([10^-2 10^3])
    set(gca, 'FontSize', 16)

if(print_to_pdf == 1)
    h = gcf;
    set(h, 'PaperOrientation', 'Landscape')
    set(h, 'PaperUnits', 'normalized')
    set(h, 'PaperPosition', [0 0 1 1])
    print(gcf, '-dpdf', '-painters', './T2S_vs_B0.pdf')
end