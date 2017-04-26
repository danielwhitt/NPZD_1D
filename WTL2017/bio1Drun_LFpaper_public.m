% initialize and timestep 1D NPZD model examples for:
% Whitt, Taylor and Levy (2017) Synoptic to planetary scale wind variab
% ility enhances phytoplankton biomass at ocean fronts
% J. Geophys. Res. Oceans

%% set the example flag below to change the scenarios

% clean up environment
clear all
close all
tic
global exampleflag

%% choose example to run: (all 3 are required to produce the figures)
exampleflag = 2; % full w (1), mean w (2), perterbation w (3)

%% specify filename where initial biogeochemical profiles are stored
inifilename = 'S020matNB2I.mat'; 

%% logic to set the number of experiments to be run (which depends on the chosen example)
totexpts = 25;
exptnums = 1:25;

%% logic to set some global variable values based on the chosen example
global sinking varW varMix nexpt
sinking = 1; % sinking detritus
varMix = 0; % constant background vertical diffusivity (0) or spatially and temporally and spatially variable vertical diffusivity (1)
varW = 1;   % timevariable vertical velocity (1) or zero vertical velocity (0)

%% define other global parameters
global dz nlev zarray wd delta zetad zetad2 sigmad sigmad2 Rm capdelta gamman alpha Vm kN kz kp I0 kback

%% nudging/restoring for deep N values
global nudgingtimescale nudgingdeepNvalue
nudgingtimescale = 0;
nudgingdeepNvalue = 5;

%% set global grid paramters
nlev = 600;
dz = .5;
z = .25:.5:299.75;
zarray = z;

%% mixing parameters; read in vertical diffusivity profiles from a file if needed
kback = 2e-5; % background vertical diffusivity


%% vertical advection parameters; read in vertical velocity time series from a file if necessary
% define vertical velocity averaging approach here:
global warray

if varW == 1
    if exampleflag == 1
        load w_filt.mat warray
    elseif exampleflag == 2
        zarray1 = zarray;
        clear zarray
        load w_filt.mat zarray oceantime
        warray = repmat(zarray(1745,:)./oceantime(1745), [1800 1]); % define mean velocity over 72.72 days
        zarray = zarray1;
        clear oceantime
    elseif exampleflag == 3
        zarray1 = zarray;
        clear zarray
        load w_filt.mat zarray warray oceantime
        warray = warray - repmat(zarray(1745,:)./oceantime(1745), [1800 1]); % subtract mean velocity over 72.72 days
        zarray = zarray1;
        clear oceantime
    end
end

%% parameters for ode45 timestepping solver
tols = 1e-6.*ones(nlev.*4,1);
tols(2*nlev+1:3*nlev) = 1e-8;
options = odeset('RelTol',1e-6,'AbsTol',tols);

%% additional time stepping parameters
dt = 1; % big step/output step size in days
nj = 75; % number of big steps
nt = 2001; % small steps
np = 51; % intermediate output steps

%% initialize output grids
Noutmat = zeros(4.*nlev,np,nj,totexpts);
Ooutmat = zeros(4.*nlev,nj,totexpts);

%% main loop (over multiple simulations)
for nexpt = exptnums 
%% specify or load the initial condition
 
 %S0 = zeros(4.*nlev,1);
 %S0(1:nlev) = 14;
 %S0(nlev+1:2*nlev) = 2;
 %S0(2*nlev+1:3*nlev) = 2;
 %S0(3*nlev+1:4*nlev) = 2;
load(inifilename,'S0'); 

S0 = S0';

%% biogeochemical parameters are set inside the main nexpt loop so they can be varied
% between experiments 
%% remineralization parameters
delta = 0.15; % inverse remineralization timescale 1/days
wd = 5; % detrital sinking velocity m/d
%% zooplankton mortality parameters
zetad = .03; % linear zooplankton mortality rate 1/d
zetad2 = 0.35; % quadratic zooplankton morality rate 1/(d mmol N/m^3)
%% phytoplankton mortality parameters
sigmad = 0.03; %linear phytoplankton senescence rate 1/d 
sigmad2 = 0;  % quadratic phytoplankton mortality rate 1/(d mmol N/m^3)
%% Grazing parameters
Rm = .5; % maximum zooplankton grazing rate 1/d
capdelta = 1.5; % Ivlev parameter 1/(mmol N/m^3)
gamman = 0.3; % excretion efficiency/sloppy feeding parameter
%% phytoplankton growth rate parameters
% light/nutrient uptake parameters:
I0=158.075; % W/m^2 surface PAR
kp = 0.02; %  self shading coefficient 1/(mmol N/m^2)
kz = 0.05; % attenuation of light due to seawater 1/m
alpha = .15; % Initial slope of P-I curve 1/(d W/m^2)
Vm = 1.0; % maximum nutrient uptake rate 1/d
kN = .1; %mmol-N/m^3 uptake half saturation constant

%% initialize time and biogeochemical arrays for timestepping
S0mat = zeros(length(S0),np,nj);
tmat = zeros(nj,1);

%% main timestepping loop
for ij = 1:nj
tmat(ij) = dt.*(ij-1);
timepts = linspace(dt.*(ij-1),dt.*ij,nt)';
[tout,Sout] = ode45(@(tout,Sout) bioodeWTL2017(tout,Sout),timepts,S0,options);
S0 = Sout(end,:);
S0mat(:,:,ij) = Sout(linspace(1,nt,np),:)';  
display(num2str(ij))
display('time:')
dt*ij
% plot the output occasionally:
if mod(dt*ij,10) ==0
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

Noutmat(:,:,:,nexpt) = squeeze(S0mat(:,:,:)); % np saves per unit dt
Ooutmat(:,:,nexpt) = squeeze(nanmean(S0mat(:,:,:),2)); % 1 save per unit dt (time-mean over np substeps)
end
pause(2)

%% Save output
display('SAVING OUTPUT...')

if exampleflag == 1
    savename = '1DNPZD_fullw';
elseif exampleflag == 2
    savename = '1DNPZD_meanw'
elseif exampleflag == 3
        savename = '1DNPZD_pertw'
end
save(strcat(savename,'_out.mat'),'Ooutmat');
save(strcat(savename,'_params.mat'),'dz','dt','tmat','zarray','delta','wd','zetad','zetad2','sigmad','sigmad2','Rm','capdelta','gamman',...
   'I0','kp','kz','alpha','Vm','kN','kback','warray');

