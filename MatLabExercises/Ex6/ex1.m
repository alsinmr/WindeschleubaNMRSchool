
clear; clc; close all

%% define constants

mu0 = 4*pi*1e-7;                    %V s / A / m
mu04pi = 1e-7;
hbar = 1.05457266e-34;              %J s / rad
h = 2*pi*hbar;                      %J s
gammaH = 26.7522128e7;              %rad / s / T
gammaN = -2.7116e7;                 %rad / s / T
rNH = 1.1*1e-10;                    %m 
omegaS = 14.1*gammaN;
omegaI = 14.1*gammaH;
deltaHN = -2*mu04pi*gammaN*gammaH*hbar/rNH^3;
sigmazz00=0;
tauc=1e-10;

%% Calculate Solution CSA=0

T1S = T1(deltaHN,sigmazz00,omegaI,omegaS,tauc);
T1I = T1(deltaHN,sigmazz00,omegaS,omegaI,tauc);
T1IS = T1x(deltaHN,omegaS,omegaI,tauc);
R = [1/T1I,1/T1IS;1/T1IS,1/T1S];
t = 0:0.01:15;
Sz=zeros(size(t));
Iz=zeros(size(t));
for k=1:length(t)  
    z = expm(-R*t(k))*[0;-2]+[10;1];
    Iz(k)=z(1);
    Sz(k)=z(2);
end

figure;
plot(t,Sz,t,exp(-t/T1S)*(-2)+1)
ylabel('<S_z(t)>')
xlabel('t [s]')
set(gca,'FontSize',14);
legend('coupled system','T_1 decay')
