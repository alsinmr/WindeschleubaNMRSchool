
clear; clc;

%% define constants

mu0 = 4*pi*1e-7;                    % V s / A / m
mu04pi = 1e-7;
hbar = 1.05457266e-34;              % J s / rad
h = 2*pi*hbar;                      % J s
gammaH = 26.7522128e7;              % rad / s / T
gammaN = -2.7116e7;                 % rad / s / T
rNH = 1.1*1e-10;                    % m 
deltaHN = -2*mu04pi*gammaN*gammaH*hbar/rNH^3;
omegaH = 14.1*gammaH;
omegaN = 14.1*gammaN;
tauc=1e-10;


%% Calculate Solution with CSA and cross-correlated cross-relaxation

sigmazzH = 100e-6;
sigmazzN = 100e-6;

T1SS = T1(deltaHN,sigmazzN,omegaH,omegaN,tauc);
T1II = T1(deltaHN,sigmazzH,omegaN,omegaH,tauc);
T1IS = T1x(deltaHN,omegaN,omegaH,tauc);
T1ISIS = T1z(deltaHN, sigmazzH, sigmazzN, omegaH, omegaN, tauc);
T1ISI = T1y(deltaHN, sigmazzH, omegaH, tauc);
T1ISS = T1y(deltaHN, sigmazzN, omegaN, tauc);

R = [1/T1II ,  1/T1IS  , 1/T1ISI ;...
     1/T1IS ,  1/T1SS  , 1/T1ISS  ;...
     1/T1ISI, 1/T1ISS , 1/T1ISIS];
 
t = 0:0.01:15;
Sz=zeros(size(t));
Iz=zeros(size(t));
IzSz = zeros(size(t));
for k=1:length(t)  
    z = expm(-R*t(k))*[0;-2;0]-[gammaH/gammaN;-1;0];
    Iz(k)=z(1);
    Sz(k)=z(2);
    IzSz(k)=z(3);
end

plot(t,exp(-t/T1SS)*(-2)+1,t,Sz)
ylabel('<S_z(t)>')
xlabel('t [s]')
set(gca,'FontSize',14);
legend('T_1 decay','coupled system')

%% Calculate Solution with CSA and without cross-correlated cross-relaxation
sigmazz00 = 100e-6;

T1SS = T1(deltaHN,sigmazz00,omegaH,omegaN,tauc);
T1II = T1(deltaHN,sigmazz00,omegaN,omegaH,tauc);
T1IS = T1x(deltaHN,omegaN,omegaH,tauc);
R = [1/T1II,1/T1IS;1/T1IS,1/T1SS];
t = 0:0.01:15;
Sz2=zeros(size(t));
Iz=zeros(size(t));
for k=1:length(t)  
    z = expm(-R*t(k))*[0;-2]+[10;1];
    Iz(k)=z(1);
    Sz2(k)=z(2);
end
figure
plot(t,exp(-t/T1SS)*(-2)+1,t,Sz2)
ylabel('<S_z(t)>')
xlabel('t [s]')
set(gca,'FontSize',14);
legend('T_1 decay','coupled system')


