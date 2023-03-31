%% simple one-spin pulse acquire simulation in Liouville space, but without relaxation
% Note: everything in SI units for clarity, frequncies in linear
% frequencies!


%% 1) clear workspace and close all figures
clear, close all

%% 2) build spin operator matrices
Ix = 1/2*[0 1; 1 0];
Iy = 1/2*[0 -1i; 1i 0];
Iz = 1/2*[1 0; 0 -1];
Ip = Ix+1i*Iy;
Im = Ix-1i*Iy;
E  = eye(2);


%% 3) set Hamiltonian parameters and build the matrix
offset = 4e3; % 20kHz offset
H0 = offset*Iz;


%% 4) set simulation parameters
rho0 = Iz; % starting density operator
pulseOp = Iy; % pulse operator
detOp = Ix + 1i*Iy; % detection operator
% detOp = Iz; % detection operator

rho_eq = Iz;

dt = 10e-6; % sampling step, "dwell time"
nPoints = 1024; % length of acquisition
t = (0:dt:(nPoints-1)*dt); % time vector

sig = zeros(1,nPoints); % pre-allocation of signal vector

%% 5) acutal simulation

%convert everything to Liouville space
H0_L = kron(H0,E)-kron(E,transpose(H0));
rho0_L = reshape(transpose(rho0),[],1);
pulseOp_L = kron(pulseOp,E)-kron(E,transpose(pulseOp));
detOp_L = reshape(transpose(detOp),[],1);
rho_eq_L = reshape(transpose(rho_eq),[],1);


% start sim
rho_L = rho0_L;

% apply pulse propagator to initial density operator
Upulse_L = expm(-1i*pi/2*pulseOp_L);
rho_L = Upulse_L*rho_L;

% build propagator of free evolution, and acquire
U0_L = expm(-1i*2*pi*H0_L*dt);

for it = 1:nPoints
   sig(it)= rho_L'*detOp_L; % detect
   rho_L = U0_L*(rho_L); % propagation
end



%% 6) plot real and imaginary part, and relaxation on top
figure(1)
title('raw signal')
clf
hold on
plot(t*1e3,real(sig),'b');
plot(t*1e3,imag(sig),'r');


xlabel('t / ms')
ylabel('signal / a.u.')
box on
legend('Re','Im','T1/T2')



%% 7) apodize, i.e. multiply with window function
twin = 1e-3; %apodization parameter
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

%% 8) Fourier transform with Zero-filling, construction of frequency vector

%fft with zerofilling two twice the original points
sig_apo(1)=sig_apo(1)/2;
spec = fftshift(fft(sig_apo,2*nPoints))*dt;

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
plot(freq/1e3,imag(spec),'r');
plot(freq/1e3,abs(spec),'--k');
xlabel('freq / kHz')
ylabel('Spectrum / a.u.')
box on
legend('Re','Im','Abs')
