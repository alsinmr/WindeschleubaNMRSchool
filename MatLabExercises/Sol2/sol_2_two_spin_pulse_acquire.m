%% simple one-spin pulse acquire simulation and signal processing
% Note: everything in SI units for clarity, frequncies in linear
% frequencies!


%% 1) clear workspace and close all figures
clear, close all

%% 2) build spin operator matrices
Ix = 1/2*[0 1; 1 0];
Iy = 1/2*[0 -1i; 1i 0];
Iz = 1/2*[1 0; 0 -1];
E  = eye(2);

%construct spin operators in two-spin basis
I1x = kron(Ix,E);
I2x = kron(E,Ix);
I1y = kron(Iy,E);
I2y = kron(E,Iy);
I1z = kron(Iz,E);
I2z = kron(E,Iz);

%% 3) set Hamiltonian parameters and build the matrix
offset1 = 4e3; % offset 1
offset2 = -3e3; % offset 2
J = 200;


H0 = offset1*I1z + offset2*I2z + J*(I1x*I2x + I1y*I2y + I1z*I2z);


%% 4) set simulation parameters
rho0 = I1z+I2z; % starting density operator
pulseOp = (I1y + I2y); % pulse operator
detOp = (I1x+I2x)+1i*(I1y+I2y) ; % detection operator

dt = 50e-6; % sampling step, "dwell time"
nPoints = 4096; % length of acquisition
t = (0:dt:(nPoints-1)*dt); % time vector

sig = zeros(1,nPoints); % pre-allocation of signal vector

%% 5) acutal simulation

rho=rho0;

% apply pulse propagator to initial density operator
Upulse = expm(-1i*pi/2*pulseOp);
rho = Upulse*rho*Upulse';

% build propagator of free evolution, and acquire
U0 = expm(-1i*2*pi*H0*dt);
for it = 1:nPoints
   sig(it)= trace(detOp*rho); % detect
   rho = U0*rho*U0'; % propagation
end

%% plot real and imaginary part
figure(1)
title('raw signal')
clf
hold on
plot(t*1e3,real(sig),'b');
plot(t*1e3,imag(sig),'r');
xlabel('t / ms')
ylabel('signal / a.u.')
box on
legend('Re','Im')

%% apodize, i.e. multiply with window function
twin = 30e-3; %apodization parameter
win = exp(-t/twin);

%apply window to signal
sig_apo = sig.*win;

figure(2)
title('apodized signal')
clf
hold on
plot(t*1e3,real(sig_apo),'b');
plot(t*1e3,imag(sig_apo),'r');
xlabel('t / ms')
ylabel('signal / a.u.')
box on
legend('Re','Im')

%% Fourier transform with Zero-filling, construction of frequency vector

%fft with zerofilling two twice the original points
sig_apo(1)=sig_apo(1)/2;
spec = fftshift(fft(sig_apo,2*nPoints));

%frequency vector
N=numel(spec);
nyqFreq = 1/(2*dt);
unitAxis = 2/N * ((0:N-1)-fix(N/2));
freq = nyqFreq * unitAxis;

figure(3)
title('spectrum')
clf
hold on
plot(freq/1e3,real(spec),'b');
xlabel('freq / kHz')
ylabel('Spectrum / a.u.')
box on
legend('Re','Im','Abs')
