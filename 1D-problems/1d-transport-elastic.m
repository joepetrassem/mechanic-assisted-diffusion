clear all
clf
close all
clc

L =1;
nbin = 50;
dx = L/(nbin);
mu0 = 0;
R = 8.31;
T = 1;%298;
c_m = 2;
D = 0.1;
V_c = 1;%6*10e-4;
sig_si = 20;
lambda = 2;

params=[L nbin dx mu0 R T c_m D V_c lambda];
M = zeros(2*nbin+1);
M(1:nbin-2,1:nbin-2) = eye(nbin-2);

options = odeset(Mass=M);

c0 = ones(nbin,1);
u0 = zeros(nbin+1,1);
w0 = [c0;u0];
% c0(1) = 0.3;
% c0(end) = 1;

[t,w]=ode15s(@(t,w) DcDt(t,w,params), [0 1.05], w0,options);

c = w(:,1:nbin);
x=linspace(dx/2,L-dx/2,nbin);

for i=1:length(t)
    figure(1); hold on; plot(x,c(i,:),'k')
end
xlabel('x')
ylabel('c')


function dcdt = DcDt(t,w,params)


L=params(1); 
nx=params(2); 
dx=params(3); 
mu0=params(4); 
R=params(5); 
T=params(6); 
c_m=params(7);
D=params(8); 
V_c=params(9); 
lambda=params(10);
sig_bound = 0;
c=w(1:nx);
u = w(nx+1:end);
dcdt = zeros(2*nx + 1,1);


x=linspace(dx/2,L-dx/2,nx);

epsilon = (u(2:end) - u(1:end-1))/dx - c.*V_c;

sigma = epsilon.*lambda;
sigma(end+1) = -sigma(end);  %bc ghost bin

dcdt(nx+1) = u(1);
for i = 1:nx  %check this,
    dcdt(nx+1+i) = sigma(i+1) - sigma(i);
end



mu = R*T*log(c./(c_m - c)) - 100*sigma(1:end-1).*x';



c_edge = (c(1:end-1) + c(2:end))/2;


fluxes = zeros(nx+1,1);
for i = 2:length(fluxes)-1
    fluxes(i) = -D/(R*T)* c_edge(i-1)*(c_m-c_edge(i-1))/c_m * (mu(i)-mu(i-1))/dx;
end
fluxes(1) = 0;
fluxes(end) = (1)*(0.5*(1+tanh(1000*(t-1))));   %turns on current at t=1


div_flux = zeros(nx,1);
for i = 1:length(div_flux)
    div_flux(i) = (fluxes(i+1) - fluxes(i))/dx;
end

for i = 1:nx
    dcdt(i) = -div_flux(i);
end



end
