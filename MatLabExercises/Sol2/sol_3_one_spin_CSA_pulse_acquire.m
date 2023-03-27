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

%% 3) set simulation parameters
rho0 = Ix; % starting density operator
detOp = Ix + 1i*Iy; % detection operator

dt = 10e-6; % sampling step, "dwell time"
nPoints = 1024; % length of acquisition
t = (0:dt:(nPoints-1)*dt); % time vector

sig = zeros(1,nPoints); % pre-allocation of signal vector

%% 4) give powder parameters
% CSA tensor, in kHz
shift_xx = -10e3;
shift_yy = -5e3;
shift_zz = 10e3;

CSA_PAS = diag([shift_xx shift_yy shift_zz]);

% load powder angles
% load leb_2ang_rank_11 % loads angles and corresponding weights
load leb_2ang_rank_101.mat
% load rep_2ang_6400pts_sph


%% 5) loop over powder angles

for i_orient = 1:numel(alphas)
    
    % rotate CSA tensor
    R = erot(alphas(i_orient),betas(i_orient),gammas(i_orient));
    CSA_LAB = R*CSA_PAS*R';
    
       
    % set Hamiltonian parameters and build the matrix
    offset = CSA_LAB(3,3); 
    H0 = offset*Iz;
    
    % build propagator of free evolution, and acquire
    U0 = expm(-1i*2*pi*H0*dt);
    %reset density matrix
    rho = rho0;
    for it = 1:nPoints
        sig(it)=sig(it)+ weights(i_orient)*trace(detOp*rho); % detect
        rho = U0*rho*U0'; % propagation
    end
    
     
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
twin = 0.2e-3; %apodization parameter
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

% Fourier transform with Zero-filling, construction of frequency vector

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
% plot(freq/1e3,imag(spec),'r');
% plot(freq/1e3,abs(spec),'--k');
xlabel('freq / kHz')
ylabel('Spectrum / a.u.')
box on
legend('Re','Im','Abs')
