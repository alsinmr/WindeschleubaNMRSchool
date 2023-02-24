%% heteronuclear two-spin pulse acquire simulation in Liouville space, with random field relaxation
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
Ia = [1 0; 0 0];
Ib = [0 0; 0 1];
E  = eye(2);

%construct spin operators in two-spin basis
I1x = kron(Ix,E);
I2x = kron(E,Ix);
I1y = kron(Iy,E);
I2y = kron(E,Iy);
I1z = kron(Iz,E);
I2z = kron(E,Iz);
I1p = kron(Ip,E);
I2p = kron(E,Ip);
I1m = kron(Im,E);
I2m = kron(E,Im);
I1a = kron(Ia,E);
I2a = kron(E,Ia);
I1b = kron(Ib,E);
I2b = kron(E,Ib);
E = eye(4);


%% 3) set Hamiltonian parameters and build the matrix
offset1 = 8e3; % offset 1
offset2 = -10e3; % offset 2
J = 400;

H0 = offset1*I1z + offset2*I2z + J*(I1z*I2z);

k_z_1 = 0.1e3;
k_xy_1 = 0.1e3;

k_xy_2 = 0.1*J;
k_z_2 = k_xy_2;

%% 4) set simulation parameters
rho0 = I1z; % starting density operator
pulseOp = I1y; % pulse operator
detOp = I1x + 1i*I1y; % detection operator
% detOp = Iz; % detection operator

rho_eq = I1z;

dt = 10e-6; % sampling step, "dwell time"
nPoints = 1024*4; % length of acquisition
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
I1x_doub_comm = ;
I1y_doub_comm = ;
I1z_doub_comm = ;

I2x_doub_comm = ;
I2y_doub_comm = ;
I2z_doub_comm = ;

%generate Relaxation superoperator
R = k_z_1*(I1z_doub_comm)+k_xy_1*(I1x_doub_comm+I1y_doub_comm)...
   +k_z_2*(I2z_doub_comm)+k_xy_2*(I2x_doub_comm+I2y_doub_comm);


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
twin = Inf; %apodization parameter
win = exp(-t/twin);

%apply window to signal
sig_apo = sig.*win;

%% 8) Fourier transform with Zero-filling, construction of frequency vector

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
plot(freq/1e3,imag(spec),'r');
plot(freq/1e3,abs(spec),'--k');
xlabel('freq / kHz')
ylabel('Spectrum / a.u.')
box on
legend('Re','Im','Abs')

xlim(offset1/1e3+[-1 1])


%% 9) Basis Transformation

%construct the full spherical Basis
clear B
B(:,1)=reshape(transpose(E),[],1);            
B(:,2)=reshape(transpose(I1z),[],1);            
B(:,3)=reshape(transpose(I2z),[],1);            
B(:,4)=reshape(transpose(I1z*I2z),[],1);      
B(:,5)=reshape(transpose(I1x),[],1);        
B(:,6)=reshape(transpose(I2x),[],1);
B(:,7)=reshape(transpose(I1y),[],1);     
B(:,8)=reshape(transpose(I2y),[],1); 
B(:,9)=reshape(transpose(I1x*I2z),[],1);     
B(:,10)=reshape(transpose(I1z*I2x),[],1);
B(:,11)=reshape(transpose(I1y*I2z),[],1);
B(:,12)=reshape(transpose(I1z*I2y),[],1);
B(:,13)=reshape(transpose(I1x*I2x),[],1);
B(:,14)=reshape(transpose(I1y*I2y),[],1);
B(:,15)=reshape(transpose(I1x*I2y),[],1);
B(:,16)=reshape(transpose(I1y*I2x),[],1);

% normalize Eigenvectors
for ii=1:size(B,2)
    B(:,ii)=B(:,ii)/norm(B(:,ii));
end

%check that B is unitary
B*B';

%convert to spherical tensor basis
H0_L_B = B'*H0_L*B;
rho0_L_B = B'*rho0_L;
pulseOp_L_B = B'*pulseOp_L*B;
detOp_L_B = B'*detOp_L;
rho_eq_L_B = B'*rho_eq_L;
R_B=B'*R*B;