clc
clear all
close all
% This code solves the nondimensional version of two-layer QG model of Phillips (1954) using the psudospectral method. Pi terms for the
%control simulation is defined as: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Backingham Pi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pi terms      Symbole          Dimensional value            Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pi1          tau_f=15           U/(lambda*aT)             %Rayleigh  friction time scale      %
% Pi2          tau_d=100          U/(lambda*aM)             %Newtonian relaxation time scale    %
% Pi3          beta=0.196         (B*lambda^2)/U            %Gradient of f(df/dy)               %
% Pi4          sigma=3.5          w/lambda                  %Width of the jet                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau_f = 15.;
tau_d = 100.;
nu    = 0.01;
beta  = 0.196;
sigma = 3.5;
p = 8.0;
%Domain size(channel)
Lx = 46;
Ly = 23;
%time step
N  = 100000;    %Total number of time step
M  = 40;       %Number of time step which is equal to 1 model day(5 model days are equal to 1 real day)
dt = 1/40;     %Time step

%Number of grid points in x and y direction
nx = 128;
ny = 192;

%%
nux = 8.0*((Lx/nx/pi)^p)/dt/8.0;
nux
nuy = 8.0*((1./ny)^p)/dt/8.0;
nuy
x = Lx*(0:nx-1)/nx;
kx = [0:nx/2 -nx/2+1:-1]*2*pi/Lx;
[D,y]=cheb(ny);     %Cheb function determines the grid points and d/dy matrix in y direction 
y=y*Ly;
D=D/Ly;
D2=D*D;
I=eye(ny+1,ny+1);
PSIR=zeros(2,nx,ny+1);
psiR=zeros(2,nx,ny+1);
%Initialization
for j=1:ny+1
    PSIR(:,:,j) = -sigma*tanh((y(j))/sigma);
    psiR(1,:,j)=fft(squeeze(PSIR(1,:,j)));
    psiR(2,:,j)=fft(squeeze(PSIR(2,:,j)));
end
psiR(:,nx/2+1,:)=0.0;

m=0;
PSIO=zeros(2,nx,ny+1);
psi=zeros(2,nx,ny+1);
%Add perturbation
for i=1:nx
    for j=1:ny+1
        PSIO(1,i,j) = 1e-2*sin(5*2*pi*x(i)/Lx)*exp(-y(j)^2/(2*sigma^2))+1e-8*randn*exp(-y(j)^2/(1.*sigma^2));
        PSIO(2,i,j) = 1e-3*sin(4*2*pi*x(i)/Lx)*exp(-y(j)^2/(2*sigma^2))+1e-9*randn*exp(-y(j)^2/(1.*sigma^2));
    end
end

PSIO(1,:,:) = squeeze(PSIR(1,:,:)+PSIO(1,:,:));

for j=1:ny+1
    psi(1,:,j) = fft(squeeze(PSIO(1,:,j)));
    psi(2,:,j) = fft(squeeze(PSIO(2,:,j)));
end
psi(:,nx/2+1,:)=0.0;
qpvo=zeros(2,nx,ny+1);
%Potential vorticity(PV) in the upper level:
%PV1=d2(PSI1)/dx2+d2(PSI1)/dy2+(-1)^1(PSI1-PSI2)+B*y=PQV1+B*y
%Potential vorticity(PV) in the lower level:
%PV2=d2(PSI2)/dx2+d2(PSI2)/dy2+(-1)^2(PSI1-PSI2)+B*y=PQV2+B*y
%So, 
%PQV1=d2(PSI1)/dx2+d2(PSI1)/dy2+(-1)^1(PSI1-PSI2) 
%PQV2=d2(PSI2)/dx2+d2(PSI2)/dy2+(-1)^2(PSI1-PSI2)
 
%PSI1 and PSI2 can be calculated from PQV1 and PQV2 equations (as it works in PSICalc.m function)
%A*(PSI1+PSI2)=(PQV1+PQV2)
%A*(PSI1-PSI2)=(PQV1-PQV2)
%where A=(d2/dx2+d2/dy2)

for i=1:nx
    qpvo(1,i,:)=(D2-(kx(i))^2*I)*squeeze(psi(1,i,:))-squeeze(psi(1,i,:)-psi(2,i,:));
    qpvo(2,i,:)=(D2-(kx(i))^2*I)*squeeze(psi(2,i,:))+squeeze(psi(1,i,:)-psi(2,i,:));
end
qpvo(:,nx/2+1,:)=0.0;
%%
dx_psi=zeros(2,nx,ny+1);
dx_qpv=zeros(2,nx,ny+1);
dy_psi=zeros(2,nx,ny+1);
dy_qpv=zeros(2,nx,ny+1);
nn=zeros(2,nx,ny+1);
QPV=zeros(2,nx,ny+1,N/M,'single');
PSI=zeros(2,nx,ny+1,N/M,'single');
dx_PSI=zeros(2,nx,ny+1);
dy_PSI=zeros(2,nx,ny+1);
dx_QPV=zeros(2,nx,ny+1);
dy_QPV=zeros(2,nx,ny+1);

tic
for n=1:N
    psi = PSICalc(qpvo,psiR,D2,kx);      %PSICalc fuction calculates streamfunction using vorticity 
    
    for i=1:nx
        dx_psi(:,i,:)=1i*kx(i)*psi(:,i,:);
        dx_qpv(:,i,:)=1i*kx(i)*qpvo(:,i,:);
        dy_psi(1,i,:)=D*squeeze(psi(1,i,:));
        dy_psi(2,i,:)=D*squeeze(psi(2,i,:));
        dy_qpv(1,i,:)=D*squeeze(qpvo(1,i,:));
        dy_qpv(2,i,:)=D*squeeze(qpvo(2,i,:));
    end
    dx_psi(:,nx/2+1,:)=0.0;
    dy_psi(:,nx/2+1,:)=0.0;
    dx_qpv(:,nx/2+1,:)=0.0;
    dy_qpv(:,nx/2+1,:)=0.0;
    
    for j=1:ny+1
        dx_PSI(1,:,j) = real(ifft(squeeze(dx_psi(1,:,j))));
        dx_PSI(2,:,j) = real(ifft(squeeze(dx_psi(2,:,j))));
        dx_QPV(1,:,j) = real(ifft(squeeze(dx_qpv(1,:,j))));
        dx_QPV(2,:,j) = real(ifft(squeeze(dx_qpv(2,:,j))));
        dy_PSI(1,:,j) = real(ifft(squeeze(dy_psi(1,:,j))));
        dy_PSI(2,:,j) = real(ifft(squeeze(dy_psi(2,:,j))));
        dy_QPV(1,:,j) = real(ifft(squeeze(dy_qpv(1,:,j))));
        dy_QPV(2,:,j) = real(ifft(squeeze(dy_qpv(2,:,j))));
    end
    %Jacobian of PV and Stream function
    %[J(Phi,PSI)=dPsi/dx*dPhi/dy-dPsi/dy*dPhi/dx)]
    %PV=QPV+beta*y
    NN = (dx_PSI.*(dy_QPV+beta)-dy_PSI.*dx_QPV);
    for j=1:ny+1
        nn(1,:,j)=fft(squeeze(NN(1,:,j)));
        nn(2,:,j)=fft(squeeze(NN(2,:,j)));
    end
    nn(:,nx/2+1,:)=0.0;
    Anew=-nn;
    for i=1:nx
        Anew(1,i,:)=squeeze(Anew(1,i,:)+(psi(1,i,:)-psi(2,i,:)-psiR(1,i,:))/tau_d);
        Anew(2,i,:)=squeeze(Anew(2,i,:)-(psi(1,i,:)-psi(2,i,:)-psiR(2,i,:))/tau_d)-(D2-(kx(i))^2*I)*squeeze(psi(2,i,:))/tau_f;
    end
    Anew(:,nx/2+1,:)=0.0;
    %Time-stepping: Adams–Bashforth methods
    if(n==1)
        qpvn=qpvo+dt*Anew;
    else
        qpvn=qpvo+0.5*dt*(3.*Anew-Aold);
    end
    Aold=Anew;
    qpvn(:,nx/2+1,:)=0.0;
    %Hyperviscosity: Ensure numerical stability and absorbs enstrophy at small scales.
    qpvo = HyperVis(qpvn,kx,nux,nuy,p,dt,n);
  
    if(mod(n,M)==0)   %save the data every model day
        
        m=m+1;
        disp(m)
        psi = PSICalc(qpvo,psiR,D2,kx);
        for j=1:ny+1
            QPV(1,:,j,m) = real(ifft(squeeze(qpvo(1,:,j))));
            QPV(2,:,j,m) = real(ifft(squeeze(qpvo(2,:,j))));
            PSI(1,:,j,m) = real(ifft(squeeze(psi(1,:,j))));
            PSI(2,:,j,m) = real(ifft(squeeze(psi(2,:,j))));
        end
        CFL = dt*max(max(max(abs(dy_PSI(:,:,:)))))/(Lx/(nx))
        toc
        tic
    end
    
end
 save(filename,'QPV','PSI','x','y','-v7.3')
