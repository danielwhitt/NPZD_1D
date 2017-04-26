function dS = bioodeWTL2017(t,S)
% ODE 1-D bio model from 
% Whitt, D. B., Taylor, J. R., and Lévy, M. and (2017),
% Synoptic-to-planetary scale wind variability enhances phytoplankton
% biomass at ocean fronts
% J. Geophys. Res. Oceans
%
% S = (N,P,Z,D) is the state vector
% each of N,P,Z,D have nlev vertical levels with spacing dz meters

global sinking varW varMix dz nlev wd delta zetad zetad2 gamman sigmad sigmad2 nudgingtimescale nudgingdeepNvalue

dS = zeros(4.*nlev,1);


% define vertical diffusivity profile
Kappa = (3600.*24).*kpp(t); % convert to m^2/d from m^2/s

% reset negative values
mask1  = S<0;
S(mask1) = 3e-5;


%N eqns:
Ix = 1;
% no flux on top and bottom:
dS(Ix) = delta.*S(3.*nlev+Ix) + gamman.*Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    -Ufn(S,Ix,dz,nlev,t).*S(nlev+Ix) ...
    + Kappa(Ix).*(-S(Ix) + S(Ix+1))./(dz.^2);
Ix = nlev;
dS(Ix) = delta.*S(3.*nlev+Ix) + gamman.*Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    -Ufn(S,Ix,dz,nlev,t).*S(nlev+Ix) ...
    + Kappa(Ix).*(S(Ix-1) - S(Ix))./(dz.^2) ...
    - nudgingtimescale.*(S(Ix)-nudgingdeepNvalue); % nudging N
Ix = (2:nlev-1)';
if varMix == 0
dS(Ix) = delta.*S(3.*nlev+Ix) + gamman.*Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    -Ufn(S,Ix,dz,nlev,t).*S(nlev+Ix)  ...
    + Kappa(Ix).*(S(Ix-1) - 2.*S(Ix) + S(Ix+1))./(dz.^2);
else
    diffflx = zeros(length(Kappa),1);
    diffflx(Ix) = .5.*(Kappa(Ix-1)+Kappa(Ix)).*(S(Ix-1) - S(Ix))./(dz);
    diffflx(end)= .5.*(Kappa(nlev-1)+Kappa(nlev)).*(S(nlev-1)-S(nlev))./dz;
    dS(Ix) = delta.*S(3.*nlev+Ix) + gamman.*Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    -Ufn(S,Ix,dz,nlev,t).*S(nlev+Ix)  ...
    + (diffflx(Ix) - diffflx(Ix+1))./(dz);
end
if varW == 1
    dS(Ix) = dS(Ix) - (wfn(Ix-1,t).*S(Ix-1) - wfn(Ix+1,t).*S(Ix+1))./(2.*dz);
end

% P eqns:
Ix = 1;
dS(nlev+Ix) = Ufn(S,Ix,dz,nlev,t).*S(nlev+Ix) - Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    -sigmad.*S(nlev+Ix) - sigmad2.*S(nlev+Ix).^2 ...
    + Kappa(Ix).*(-S(nlev+Ix) + S(nlev+Ix+1))./(dz.^2);
Ix = nlev;
dS(nlev+Ix) = Ufn(S,Ix,dz,nlev,t).*S(nlev+Ix) - Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    -sigmad.*S(nlev+Ix) - sigmad2.*S(nlev+Ix).^2 ...
    + Kappa(Ix).*(S(nlev+Ix-1) - S(nlev+Ix))./(dz.^2);
Ix = (2:nlev-1)';
if varMix == 0
dS(nlev+Ix) = Ufn(S,Ix,dz,nlev,t).*S(nlev+Ix) - Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    -sigmad.*S(nlev+Ix) - sigmad2.*S(nlev+Ix).^2 ...
    + Kappa(Ix).*(S(nlev+Ix-1) - 2.*S(nlev+Ix) + S(nlev+Ix+1))./(dz.^2);
else
    diffflx = zeros(size(Kappa));
    diffflx(Ix) = .5.*(Kappa(Ix-1)+Kappa(Ix)).*(S(nlev+Ix-1) - S(nlev+Ix))./(dz);
    diffflx(end)= .5.*(Kappa(nlev-1)+Kappa(nlev)).*(S(nlev+nlev-1)-S(nlev+nlev))./dz;
    dS(nlev+Ix) = Ufn(S,Ix,dz,nlev,t).*S(nlev+Ix) - Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
        -sigmad.*S(nlev+Ix) - sigmad2.*S(nlev+Ix).^2 ...
    + (diffflx(Ix) - diffflx(Ix+1))./(dz);
end

if varW == 1
    dS(nlev+Ix) = dS(nlev+Ix) - (wfn(Ix-1,t).*S(nlev+Ix-1) - wfn(Ix+1,t).*S(nlev+Ix+1))./(2.*dz);
end


% Z eqns:
Ix = 1;
dS(2.*nlev+Ix) = (1 - gamman).*Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    - zetad.*S(2.*nlev+Ix) - zetad2.*S(2.*nlev+Ix).^2 ...
    + Kappa(Ix).*(-S(2.*nlev+Ix) + S(2.*nlev+Ix+1))./(dz.^2);
Ix = nlev;
dS(2.*nlev+Ix) = (1 - gamman).*Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    - zetad.*S(2.*nlev+Ix) - zetad2.*S(2.*nlev+Ix).^2 ...
    + Kappa(Ix).*(S(2.*nlev+Ix-1) - S(2.*nlev+Ix))./(dz.^2);
Ix = (2:nlev-1)';
if varMix == 0
dS(2.*nlev+Ix) = (1 - gamman).*Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
    - zetad.*S(2.*nlev+Ix) - zetad2.*S(2.*nlev+Ix).^2 ...
    + Kappa(Ix).*(S(2.*nlev+Ix-1) - 2.*S(2.*nlev+Ix) + S(2.*nlev+Ix+1))./(dz.^2);
else
    diffflx = zeros(size(Kappa));
    diffflx(Ix) = .5.*(Kappa(Ix-1)+Kappa(Ix)).*(S(2.*nlev+Ix-1) - S(2.*nlev+Ix))./(dz);
    diffflx(end)= .5.*(Kappa(nlev-1)+Kappa(nlev)).*(S(2.*nlev+nlev-1)-S(2.*nlev+nlev))./dz;
    dS(2.*nlev+Ix) = (1 - gamman).*Gfn(S(nlev+Ix)).*S(2.*nlev+Ix) ...
        - zetad.*S(2.*nlev+Ix) - zetad2.*S(2.*nlev+Ix).^2 ...
            + (diffflx(Ix) - diffflx(Ix+1))./(dz);
end
if varW == 1
    dS(2.*nlev+Ix) = dS(2.*nlev+Ix) - (wfn(Ix-1,t).*S(2.*nlev+Ix-1) - wfn(Ix+1,t).*S(2.*nlev+Ix+1))./(2.*dz);
end

% D eqns
if sinking == 1
Ix = 1;
dS(3.*nlev+Ix) = sigmad.*S(nlev+Ix) + sigmad2.*S(nlev+Ix).^2 ...
    + zetad.*S(2.*nlev+Ix) + zetad2.*S(2.*nlev+Ix).^2 ...
    - delta.*S(3.*nlev+Ix) ...
    + wd.*(0 - wdfn(Ix+1).*S(3.*nlev+Ix+1))./(dz)...
    + Kappa(Ix).*(-S(3.*nlev+Ix) + S(3.*nlev+Ix+1))./(dz.^2);
Ix = nlev;
dS(3.*nlev+Ix) = sigmad.*S(nlev+Ix) +sigmad2.*S(nlev+Ix).^2 ...
    + zetad.*S(2.*nlev+Ix) + zetad2.*S(2.*nlev+Ix).^2 ...
    - delta.*S(3.*nlev+Ix) ...
    + wd.*(wdfn(Ix-1).*S(3.*nlev+Ix-1) - wdfn(Ix).*S(3.*nlev+Ix))./(dz)...
    + Kappa(Ix).*(S(3.*nlev+Ix-1) - S(3.*nlev+Ix))./(dz.^2);
Ix = (2:nlev-1)';
if varMix == 0
dS(3.*nlev+Ix) = sigmad.*S(nlev+Ix) +sigmad2.*S(nlev+Ix).^2 ...
    + zetad.*S(2.*nlev+Ix) + zetad2.*S(2.*nlev+Ix).^2 ...
    - delta.*S(3.*nlev+Ix) ...
    + wd.*(wdfn(Ix-1).*S(3.*nlev+Ix-1) - wdfn(Ix+1).*S(3.*nlev+Ix+1))./(2.*dz)...
    + Kappa(Ix).*(S(3.*nlev+Ix-1) - 2.*S(3.*nlev+Ix) + S(3.*nlev+Ix+1))./(dz.^2);
else
    diffflx = zeros(size(Kappa));
    diffflx(Ix) = .5.*(Kappa(Ix-1)+Kappa(Ix)).*(S(3.*nlev+Ix-1) - S(3.*nlev+Ix))./(dz);
    diffflx(end)= .5.*(Kappa(nlev-1)+Kappa(nlev)).*(S(3.*nlev+nlev-1)-S(3.*nlev+nlev))./dz;
    dS(3.*nlev+Ix) = sigmad.*S(nlev+Ix) +sigmad2.*S(nlev+Ix).^2 ...
        + zetad.*S(2.*nlev+Ix) + zetad2.*S(2.*nlev+Ix).^2 ...
        - delta.*S(3.*nlev+Ix) ...
        + wd.*(wdfn(Ix-1).*S(3.*nlev+Ix-1) - wdfn(Ix+1).*S(3.*nlev+Ix+1))./(2.*dz)...
        + (diffflx(Ix) - diffflx(Ix+1))./(dz);
end
else
    Ix = 1;
dS(3.*nlev+Ix) = sigmad.*S(nlev+Ix) + sigmad2.*S(nlev+Ix).^2 ...
    + zetad.*S(2.*nlev+Ix) + zetad2.*S(2.*nlev+Ix).^2 ...
    - delta.*S(3.*nlev+Ix) ...
    + Kappa(Ix).*(-S(3.*nlev+Ix) + S(3.*nlev+Ix+1))./(dz.^2);
Ix = nlev;
dS(3.*nlev+Ix) = sigmad.*S(nlev+Ix) + sigmad2.*S(nlev+Ix).^2 ...
    + zetad.*S(2.*nlev+Ix) + zetad2.*S(2.*nlev+Ix).^2 ...
    - delta.*S(3.*nlev+Ix) ...
    + Kappa(Ix).*(S(3.*nlev+Ix-1) - S(3.*nlev+Ix))./(dz.^2);
Ix = (2:nlev-1)';
if varMix == 0
dS(3.*nlev+Ix) = sigmad.*S(nlev+Ix) + sigmad2.*S(nlev+Ix).^2 ...
    + zetad.*S(2.*nlev+Ix) + zetad2.*S(2.*nlev+Ix).^2 ...
    - delta.*S(3.*nlev+Ix) ...
    + Kappa(Ix).*(S(3.*nlev+Ix-1) - 2.*S(3.*nlev+Ix) + S(3.*nlev+Ix+1))./(dz.^2);
else
    diffflx = zeros(size(Kappa));
    diffflx(Ix) = .5.*(Kappa(Ix-1)+Kappa(Ix)).*(S(3.*nlev+Ix-1) - S(3.*nlev+Ix))./(dz);
    diffflx(end)= .5.*(Kappa(nlev-1)+Kappa(nlev)).*(S(3.*nlev+nlev-1)-S(3.*nlev+nlev))./dz;
    dS(3.*nlev+Ix) = sigmad.*S(nlev+Ix) + sigmad2.*S(nlev+Ix).^2 ...
        + zetad.*S(2.*nlev+Ix) + zetad2.*S(2.*nlev+Ix).^2 ...
        - delta.*S(3.*nlev+Ix) ...
        + (diffflx(Ix) - diffflx(Ix+1))./(dz);
end
end
if varW == 1
    dS(3.*nlev+Ix) = dS(3.*nlev+Ix) - (wfn(Ix-1,t).*S(3.*nlev+Ix-1) - wfn(Ix+1,t).*S(3.*nlev+Ix+1))./(2.*dz);
end
end

function G = Gfn(P)
% Ivlev zooplankton grazing function
global Rm capdelta 
G = Rm.*(1-exp(-capdelta.*P));
end

function U = Ufn(S,idx,dz,nlev,t)
% phytoplankton growth rate function
global alpha Vm kN 
U = ((Vm.*S(idx))./(kN+S(idx)))...
    .*(alpha.*Ifn(S,idx,dz,nlev,t)./(sqrt(Vm.^2 + (alpha.^2).*Ifn(S,idx,dz,nlev,t).^2)));
end

function I = Ifn(S,idx,dz,nlev,t)
% light function
global  kz kp I0
selfshading = 1; 
if selfshading == 1
    if idx == nlev
        I = I0.*exp(-kz.*zfn(idx) - .5.*dz.*kp.*S(nlev+1) ...
            - kp.*trapz(zfn(1:idx),S(nlev+1:nlev+idx)));
    elseif idx == 1
        I = I0.*exp(-kz.*zfn(idx) - .5.*dz.*kp.*S(nlev+1));
    else
        I = I0.*exp(-.5.*dz.*kp.*S(nlev+1)).*exp(-kz.*zfn(idx) ...
            - kp.*cumtrapz(zfn(idx),S(nlev+idx)));
    end
else
    I = I0.*exp(-kz.*zfn(idx));
end
end

function z = zfn(idx)
global zarray
z = zarray(idx)';
end


function k = kpp(t)
% vertical diffusivity profile, output in m^2/s
global zarray kback varMix Kz nexpt
if varMix == 0
k = kback.*ones(size(zarray));
k = k';
else
k = Kz(:,round(t.*10)+1,nexpt);
end
end


function w = wfn(idx,t)
% vertical velocity function (does not vary in space)
global  warray nexpt nlev varW
if varW == 1
    w = warray(max(round(t.*24),1),nexpt);
else
    w = 0;
end
% vertical structure is constant
w1 = ones(nlev,1);
w = w.*w1(idx);
end


function w = wdfn(idx)
% vertical structure function for detrital sinking 
%- here detrital sinking is constant in space and time
global nlev
w1 = ones(nlev,1);
w = w1(idx);
end
