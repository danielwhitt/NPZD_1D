% script to run an example experiment using the 1D model in 
% Whitt, D. B., Lévy, M. and Taylor, J. R. (2016), Low- and high-frequency oscillatory winds synergistically enhance nutrient entrainment and phytoplankton at fronts. 
% J. Geophys. Res. Oceans. Accepted Author Manuscript. doi:10.1002/2016JC012400

clear all;
nproc1=1 
totexpts = 1;
exptnums = (nproc1-1)*25+(1:totexpts)
%% global parameters/flags
global sinking tides varMix MLD
sinking = 1; % sinking detritus
varMix = 1; % variable vertical diffusion
MLD = 90;
tides = 0;
Nthreads = 1;
maxNumCompThreads(Nthreads)
tic

%% global grid paramters
global dz nlev zarray
nlev = 300;
dz = 1;
z=.5:dz:299.5;
zarray = z;

%% parameters for ode solve
tols = 1e-7.*ones(nlev.*4,1);
tols(2*nlev+1:3*nlev) = 1e-9;
options = odeset('RelTol',1e-7,'AbsTol',tols);
nj = 101;
nout = nj;
nt = 2001;
np = 11;
dt =1; % 1 day


%% global variables for NPZD model
global nudgingtimescale nudgingdeepNvalue
nudgingtimescale = 0;
nudgingdeepNvalue = 4;

global wd delta zetad zetad2 sigmad sigmad2 Rm capgamma gamman alpha Vm kN kz kp I0 wt omegat kback

%% output profiles
Noutmat = zeros(4.*nlev,np,nout,totexpts);
Ooutmat = zeros(4.*nlev,nout,totexpts);

for nexpt = exptnums 
%% load initial condition
 load S0_1D_schematic_HF.mat
 S0 = S0';

%% remineralization parameters
delta = .15; %1/day
wd=5; %m/day

%% zooplankton and phytoplankton mortality parameters
zetad = .03; % 1/d linear zooplankton mortality
zetad2 = 0.35; % 1/(d mmol N/m^3) quadratic zooplankton mortality

sigmad = 0.03; % 1/d phytoplankton senescence
sigmad2 = 0; % quadratic phytoplankton mortality

%% Grazing parameters
% Grazing parameters
Rm = .5; % max grazing rate
capgamma = 1.5; % Ivlev parameter
gamman = 0.3; %excretion efficiency, same as Franks 86

%% phytoplankton growth rate parameters
% light/nutrient uptake parameters:
I0 = 158.075; % W/m^2 surface irradiance % B&H:ranges from 100-300
kp = 0.02; %  m^2 /mmol-N self shading coef B&H:.03 - .07 
kz = .05; %  
% P/I slope dependence:
 %d-1 m^2/W Initial slope of P-I curve B&H: .01 - .19
alpha = .15; 
Vm = 1.0; % maximum nutrient uptake rate 1/d
kN = .1; % mmol-N/m^3 uptake half saturation %B&H: .1-.5, 


%% vertical advection parameters
wt = 0; % velocity maximum is at pi/2 (150 m) we are setting the velocity
omegat = 1.92; % omega = 8E-5

%% vertical diffusion parameters
kback = 2e-5; % vertical diffusivity m^2/s

%% initialize solution array
S0mat = zeros(length(S0),np,nj);
tmat = zeros(nj,1);

%% time step
for ij = 1:nj
tmat(ij) = dt.*(ij-1);
timepts = linspace(dt.*(ij-1),dt.*ij,nt)';
[tout,Sout] = ode45(@(tout,Sout) bioodeWLT2017(tout,Sout),timepts,S0,options);
S0 = Sout(end,:);
S0mat(:,:,ij) = Sout(linspace(1,nt,np),:)';  
display(num2str(ij))
display('time:')
dt*ij
%% plotting
if mod(dt*ij,10) ==2 
    clf; cla;
subplot(2,4,1),pcolor(tmat(1:ij),z,squeeze(S0mat(0.*nlev+1:1.*nlev,1,1:ij)));
shading flat
colorbar
title('N')
subplot(2,4,2),pcolor(tmat(1:ij),z,squeeze(S0mat(1.*nlev+1:2.*nlev,1,1:ij)));
shading flat
colorbar

title('P')
subplot(2,4,3),pcolor(tmat(1:ij),z,squeeze(S0mat(2.*nlev+1:3.*nlev,1,1:ij)));
shading flat
colorbar

title('Z')
subplot(2,4,4),pcolor(tmat(1:ij),z,squeeze(S0mat(3.*nlev+1:4.*nlev,1,1:ij)));
shading flat
colorbar

title('D')
subplot(2,4,5),pcolor(linspace(dt*(ij-1),dt*ij,np),z,squeeze(S0mat(0.*nlev+1:1.*nlev,:,ij)));
shading flat
colorbar
title('N')
subplot(2,4,6),pcolor(linspace(dt*(ij-1),dt*ij,np),z,squeeze(S0mat(1.*nlev+1:2.*nlev,:,ij)));
shading flat
colorbar

title('P')
subplot(2,4,7),pcolor(linspace(dt*(ij-1),dt*ij,np),z,squeeze(S0mat(2.*nlev+1:3.*nlev,:,ij)));
shading flat
colorbar
title('Z')
subplot(2,4,8),pcolor(linspace(dt*(ij-1),dt*ij,np),z,squeeze(S0mat(3.*nlev+1:4.*nlev,:,ij)));
shading flat
colorbar

title('D')
pause(.1)
end

toc
   
end
display('done:')
display(num2str(nexpt))

Noutmat(:,:,:,nexpt) = squeeze(S0mat(:,:,:));
Ooutmat(:,:,nexpt) = squeeze(nanmean(S0mat(:,:,:),2));
end
