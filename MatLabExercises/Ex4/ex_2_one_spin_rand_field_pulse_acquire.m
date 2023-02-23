%% simple one-spin pulse acquire simulation in Liouville space, with random field relaxation
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

k_z = 5e2;
k_xy = 1e2;


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
H0_L = ;
rho0_L = ;
pulseOp_L = ;
detOp_L = ;
rho_eq_L = ;

%generate the double commutators
Ix_doub_comm = ;
Iy_doub_comm = ;
Iz_doub_comm = ;

%generate Relaxation superoperator
R = k_z*(Iz_doub_comm)+k_xy*(Ix_doub_comm+Iy_doub_comm);



% start sim
rho_L = rho0_L;

% apply pulse propagator to initial density operator
Upulse_L = expm(-1i*pi/2*pulseOp_L);
rho_L = ;

% build propagator of free evolution, and acquire
U0_L = expm(-1i*2*pi*H0_L*dt-R*dt);

for it = 1:nPoints
   sig(it)= ; % detect
   rho_L = ; % propagation
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
legend('Re','Im')



%% 7) apodize, i.e. multiply with window function
twin = Inf; %apodization parameter, not needed here!
win = exp(-t/twin);

%apply window to signal
sig_apo = sig.*win;


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

%% 9) Basis Transformation

% basis transformation to operator Basis
B(:,1)=reshape(transpose(E),[],1);
B(:,2)=reshape(transpose(Iz),[],1);
B(:,3)=reshape(transpose(Ix),[],1);
B(:,4)=reshape(transpose(Iy),[],1);

% normalize Eigenvectors
for ii=1:size(B,1)
    B(ii,:)=B(ii,:)/norm(B(ii,:));
end

%check that B is unitary
B*B';


% convert to Operator Basis
%convert to spherical tensor basis
H0_L_B = B'*H0_L*B;
rho0_L_B = B'*rho0_L;
pulseOp_L_B = B'*pulseOp_L*B;
detOp_L_B = B'*detOp_L;
rho_eq_L_B = B'*rho_eq_L;
R_B=B'*R*B;