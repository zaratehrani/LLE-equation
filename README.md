# LLE-equation
%LLE Equation solver
%This code reproduces the results of Optics Communications 312(2014)134–136 

%July 26, 2021

close all
clear all

nF=512;  %Numnber of modes in fourier space 

Q=3e9;
Omega0=2*pi*(3*10^8)/(1560.5e-9);
kappa=Omega0/Q;
Time=kappa/2*.002; %This is the total time of the simulation time (old time Time=4000); 
dt=0.05;  %time step
M=3;  %30 This is the fast time, the amount of time simulated inside the cavity because we relate it to # of modes

iter=16;   %number of iterations for solving fast time steps
dtRK=dt/iter; %RK time step
N=round(Time/dt);  % This is the total time of the simulation
m=M/nF;  %fast time step
n=(-nF/2:1:nF/2-1)';  %mode indices for doing fourier transform, I believe these will double for the comb teeth indices
nt=n*m;  %These the time values inside the cavity
k=2*pi*n/M;  %these are the wave numbers (frequencies) for the modes of the FFT and comb teeth

%del_omega=2*pi*13.36e9;
%radius=2.5e-03;
%l=2*pi*radius;


D2=1;      %This is the second order dispersion term but it will be left as unity at all frequencies for now
theta=0*3;   %This is the Detuning in the frequency domain this should be adjusted to D1

%Now set some simulation parameters
E_in=1.8*sqrt(2);%4*sqrt(2);%1.01*sqrt(2),,1.2*sqrt(2) This is the amplitude at the central cavity mode in the frequency domain
g=1;       %This is the nonlinear coefficient for the FWM term
%Initial conditions: We will seed all the frequency modes with random
%amplitudes
u0=rand(length(nt),1)*(1e-5).*exp(1i*rand(length(nt),1)*2*pi);
%This is in the frequency domain so we will FFT to the time domain.
u0=ifft(fftshift(u0));


%Now create matrix to represent input at central frequency
fE_in=repmat(E_in,nF,1);
fE_in=fftshift(fft(fE_in)); %Power spectrum at central frequency
u = u0;

%   P=-1i*D2*k.^2-(1+1i*theta);
 P=-1i*.00625*n.^2-(1);
 expP=exp(P*dt/2);  
  %Now add the pump contribution
 PumpC=fE_in.*(expP-1)./P;
    
%Create empty matrix with all the data for every cavity round trip
U=zeros(nF,N);  %Time domain
V=zeros(nF,N);      %Frequency domain  
    
        %Here we are using the split step fourier method, this involves
    %1.Propagation in the frequeny domain for half the slow time step then
    %2. calculating the FWM contribution in the fast time 
    %3. Propagation of the other half in the frequency domain
    %ii
    %propagtion in frequency domain
tic %start calculating calculation time
v0=zeros(length(u),1);
v=zeros(length(u),1);
for ii=1:1:N %start looping in slow time
    %Now we add this to the starting cavity conditions, currently 
    v0 = fftshift(fft(u));
    v0 = PumpC + expP.*v0;
    %Now FFT to time domain to add FWM terms
    v = ifft(fftshift(v0));
    %Solve for FWM in time domain using runge kutta 2
    for ll=1:1:iter
        f1 = v + 0.5*dtRK*g*1i*abs(v).^2.*v;
        f2 = g*1i*abs(f1).^2.*f1;
        v = v +dtRK*f2;
    end
    %Now do other half of the propagation
    v = fftshift(fft(v)); %Go back to frequency domain
    v = PumpC + expP.*v;
    u = ifft(fftshift(v));
    U(:,ii)=u;
    V(:,ii)=v;
      
end

toc
% %% Time Domain surface plot
figure(1);clf
xscale = 100;
yscale = 2;
s = surf(abs(U(1:yscale:end,1:xscale:end)').^2);
s.EdgeColor = 'none';
set(gca, 'xtick', 1:nF/(10*yscale):nF)
set(gca, 'xticklabels', -M/2:M/(10):M/2)
set(gca, 'ytick', 1:N/(10*xscale):N/xscale)
set(gca, 'yticklabels', 0:Time/(10):Time)
axis tight
xlabel('\xi')
ylabel('t')
zlabel('|E|^2')
clear xscale yscale
view(3)



%% Frequency Domain surface plot
xscale=1;
figure(2);clf
s = surface(abs(V(1:2:end,1:xscale:end)').^2);
  s.EdgeColor = 'none';
set(gca, 'xtick', 1:nF/(16*2):nF)
set(gca, 'xticklabels', n(1):nF/(16):n(end))
set(gca, 'ytick', 1:N/(8*xscale):N/xscale)
set(gca, 'yticklabels', 0:Time/8:Time)
axis tight
xlabel('\mu')
ylabel('t')
zlabel('|E|^2')
view(3)

figure(4)
plot(1:N,smooth(sum(abs(V).^2,1),1000))

%  figure(3)
% hold on
% % plot(n(1):1:n(end),abs(V(:,1)).^2);
% % plot(n(1):1:n(end),abs(V(:,N/4)).^2);
% plot(n(1):1:n(end),abs(V(:,N/2)).^2);
%  plot(n(1):1:n(end),abs(V(:,3*N/4)).^2);
% plot(n(1):1:n(end),abs(V(:,end)).^2);
% hold off
% legend
Norm=0.5*2*pi*64e3; %Via the chemo paper
figure(5)
semilogy(n(1):1:n(end),abs(V(:,end)).^2/Norm);
ylim([10^-9 10]) 
xlim([-100 100]) 
xlabel('Mode Number \mu')
ylabel('|E|^2 Log Scale')