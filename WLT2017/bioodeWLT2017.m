function dS = bioodeWLT2017(t,S)
% ODE 1-D bio model from 
% Whitt, D. B., Lévy, M. and Taylor, J. R. (2016), Low- and high-frequency oscillatory winds synergistically enhance nutrient entrainment and phytoplankton at fronts. 
% J. Geophys. Res. Oceans. Accepted Author Manuscript. doi:10.1002/2016JC012400
% S = (N,P,Z,D) 
% each of N,P,Z,D have nlev vertical levels with spacing dz meters
% N eqns:
% S is state variable
% dSdt:
global sinking tides varMix dz nlev wd delta zetad zetad2 gamman sigmad sigmad2 wt omegat kback nudgingtimescale nudgingdeepNvalue

dS = zeros(4.*nlev,1);


if varMix == 1
    Kappa = (3600.*24).*kpp(t);
    %if t < 1
    %Kappa(Kappa > 1E-3.*86400) = 1E-3.*86400;
    %end
else
    Kappa = kback.*(3600.*24).*ones(nlev,1);
end

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
if tides == 1
    dS(Ix) = dS(Ix) - (wfn(Ix-1,wt,t,omegat).*S(Ix-1) - wfn(Ix+1,wt,t,omegat).*S(Ix+1))./(2.*dz);
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

if tides == 1
    dS(nlev+Ix) = dS(nlev+Ix) - (wfn(Ix-1,wt,t,omegat).*S(nlev+Ix-1) - wfn(Ix+1,wt,t,omegat).*S(nlev+Ix+1))./(2.*dz);
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
if tides == 1
    dS(2.*nlev+Ix) = dS(2.*nlev+Ix) - (wfn(Ix-1,wt,t,omegat).*S(2.*nlev+Ix-1) - wfn(Ix+1,wt,t,omegat).*S(2.*nlev+Ix+1))./(2.*dz);
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
if tides == 1
    dS(3.*nlev+Ix) = dS(3.*nlev+Ix) - (wfn(Ix-1,wt,t,omegat).*S(3.*nlev+Ix-1) - wfn(Ix+1,wt,t,omegat).*S(3.*nlev+Ix+1))./(2.*dz);
end
end

function G = Gfn(P)
global Rm capgamma 
%Mayzaud Poulet 78
%G = Rm.*capgamma.*P.*(1-exp(-capgamma.*P));
% Ivlev
G = Rm.*(1-exp(-capgamma.*P));
end

function U = Ufn(S,idx,dz,nlev,t)
global alpha Vm kN 
% Spitz et al. 2003
U = ((Vm.*S(idx))./(kN+S(idx)))...
    .*(alpha.*Ifn(S,idx,dz,nlev,t)./(sqrt(Vm.^2 + (alpha.^2).*Ifn(S,idx,dz,nlev,t).^2)));
end

function I = Ifn(S,idx,dz,nlev,t)
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
global zarray MLD kback varMix
if varMix == 0
k = kback.*ones(size(zarray));
else
if mod(t,8) < 1
k = 1e-2.*(.5-.5.*tanh((zarray-MLD)./3))+kback;
else
k = kback.*ones(size(zarray));
end
end
k = k';
end

 

function w = wfn(idx,wt,t,omegat)
% vertical advection function of all constituents
global nlev warray mixingflag

if mixingflag == 0
%w1 = wt.*sin(2.*pi.*omegat.*t).*sin(pi.*linspace(0,1,nlev))';
%w1 = (wt.*cos(2.*pi.*omegat.*t)+(6.737./.7071).*cos(2.*pi.*.1238.*t)).*sin(pi.*linspace(0,1,nlev))'; %gotta divide
%w1 = (wt.*cos(2.*pi.*omegat.*t)+(2./.7071).*cos(2.*pi.*.02.*t)).*sin(pi.*linspace(0,1,nlev))'; %gotta divide
%elseif wt == 2 % downwelling
%    w1 = -.2.*sin(pi.*linspace(0,1,nlev))';
if wt == 11 % upwelling
    w1 = (.1.*sin(pi.*linspace(0,1,nlev)).^0)';
 elseif wt == 0 % control, nothing
    w1 = .1.*sin(pi.*linspace(0,1,nlev))';
elseif wt == 1 % upwelling + LF mpo2   
   % w1 = ((-.2+(15./.866).*cos(2.*pi.*.1238.*t)).*sin(pi.*linspace(0,1,nlev)))';
      w1 = ((.1+(9).*cos((2.*pi.*.1238.*t)-pi./2)).*sin(pi.*linspace(0,1,nlev)).^0)';
elseif wt == 2 % upwelling plus LF
    w1 = ((.1+(9).*cos((2.*pi.*.1238.*t))).*sin(pi.*linspace(0,1,nlev)).^0)';
elseif wt == 3 % upwelling + LF 
    %w1 = ((.2+(120./.866).*cos(2.*pi.*1.1.*t)).*sin(pi.*linspace(0,1,nlev)))';
    w1 = (((9).*cos((2.*pi.*.1238.*t)-pi./2)).*sin(pi.*linspace(0,1,nlev)).^0)';
elseif wt == 4 % upwelling + LF
    w1 = (((9).*cos(2.*pi.*.1238.*t)).*sin(pi.*linspace(0,1,nlev)).^.0)';
 %   w1 = ((.2+(15./.866).*cos(2.*pi.*.1238.*t)).*sin(pi.*linspace(0,1,nlev)))';
elseif wt == 5 % downwelling + HF   
    %w1 = ((-.2+(120./.866).*cos(2.*pi.*1.1.*t)).*sin(pi.*linspace(0,1,nlev)))';
    w1 = (((8).*cos((2.*pi.*.1238.*t))).*sin(pi.*linspace(0,1,nlev)).^0)';
elseif wt == 7 % LF
    w1 = (((6./.866).*cos(2.*pi.*.1238.*t)).*sin(pi.*linspace(0,1,nlev)).^0)';
elseif wt == 9 % upwelling + LFHF
    w1 = ((.2+(15./.866).*cos(2.*pi.*.1238.*t)+(120./.866).*cos(2.*pi.*1.1.*t)).*sin(pi.*linspace(0,1,nlev)))';
elseif wt == 8 % downwelling + LFHF   
    w1 = ((-.2+(15./.866).*cos(2.*pi.*.1238.*t)+(120./.866).*cos(2.*pi.*1.1.*t)).*sin(pi.*linspace(0,1,nlev)))';
elseif wt == 10 % HF
    w1 = (((120./.866).*cos(2.*pi.*1.1.*t)).*sin(pi.*linspace(0,1,nlev)))';
elseif wt == 11 %LFHF
    w1 = (((15./.866).*cos(2.*pi.*.1238.*t)+(120./.866).*cos(2.*pi.*1.1.*t)).*sin(pi.*linspace(0,1,nlev)))';
end
w = w1(idx);

elseif mixingflag == 2 || mixingflag ==4 || mixingflag == 5 || mixingflag == 6
    w = warray(max(round(t.*24),1));
end

%w1 = (wt.*cos(2.*pi.*omegat.*t)).*sin(pi.*linspace(0,1,nlev))'; %gotta divide

%w1 = (wt.*cos(2.*pi.*omegat.*t)).*sin(pi.*linspace(0,1,nlev))'; %gotta divide
% by .8660 to get the RMS displacements in the right place -see powellrun.m
%w1 = wt.*sin(pi.*linspace(0,4,nlev))';
%w1(151:600) =0;
end


function w = wdfn(idx)
% vertical variability of detrital sinking
global nlev
w1 = ones(nlev,1);
%test = linspace(0,-300,nlev)./50 +75./50;
%w1 = (.5-tanh(test)./2)';
w = w1(idx);
end
