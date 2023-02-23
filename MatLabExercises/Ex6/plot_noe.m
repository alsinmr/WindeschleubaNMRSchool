clear; close all

%% Define Constants
gammaH = 26.7522128e7;              %rad / s / T
gammaN = -2.7116e7;                 %rad / s / T
omegaS = 14.1*gammaN;
omegaI = 14.1*gammaH;
mu0 = 4*pi*1e-7;                    % V s / A / m
mu04pi = 1e-7;
hbar = 1.05457266e-34;              %J s / rad
h = 2*pi*hbar;                      %J s
rNH = 1.1*1e-10;                    %m 

%% define eta 
tauc = 10.^(-12:0.01:-5);
delta = -2*mu04pi*gammaN*gammaH*hbar/rNH^3;
sigmazz00=0;

T1S = T1(delta,sigmazz00,omegaI,omegaS,tauc);
T1IS = T1x(delta,omegaS,omegaI,tauc);

eta = T1S./T1IS*gammaH/gammaN;

%% plot 
semilogx(omegaS*tauc,eta,'b');
ylabel('\eta')
set(gca,'FontSize',14);
axis([-1e3 -1e-3 -6 1]);
set(gca,'XDir','Reverse');
set(gca,'FontSize',14);
xlabel('\omega_0^{(S)}\tau_c')