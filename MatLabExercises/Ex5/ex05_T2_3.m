%% Problem Set 5: Relaxation Rate Constants T2
%Problem 2.3

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

B0 = ;                       %T
sigma_zz = ;                        %ppm
tauC1 = ;                      %s
tauC2 = ;                      %s
tauC3 = ;                       %s

omegaI = ;                 %rad/s, 1H Larmor frequency
omegaS = ;                 %rad/s, 15N Larmor frequency

deltaIS = -2*mu0*gammaH*gammaN*hbar/(4*pi)/rNH^3;

%% 4) Calculate T2S
%s. Eq. (4) on the exercise sheet
    %note that sigma_zz is given in ppm above
    %Hint: Elementwise multiplication in Matlab: array1.*array2
    
R2S_tauC1 = ;
R2S_tauC2 = ;
R2S_tauC3 = ;

T2S_tauC1 = 1./R2S_tauC1;
T2S_tauC2 = 1./R2S_tauC2;
T2S_tauC3 = 1./R2S_tauC3;
   
%% 5 Semi-Logarithmic Plot of T1S vs B0

figure
    semilogy(B0, T2S_tauC1)
    hold on
    semilogy(B0, T2S_tauC2)
    semilogy(B0, T2S_tauC3)
    xlabel("B_0 [T]")
    ylabel("T_{2, S} [s]")
    lgd = legend([num2str(tauC1); num2str(tauC2); num2str(tauC3)]);
        lgd.Title.String = "\tau_c";
        lgd.Location = 'Northeast';
    xlim([B0(1) B0(end)])
    ylim([10^-2 10^3])
 