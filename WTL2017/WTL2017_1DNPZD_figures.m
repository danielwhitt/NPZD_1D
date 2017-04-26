% generate some figures for: 
% Whitt, Taylor and Levy (2017) Synoptic to planetary scale wind variab
% ility enhances phytoplankton biomass at ocean fronts
% J. Geophys. Res. Oceans

close all
clear all
%% compare initial conditions, Figure S1
% ROMS, P,Z,D
load COSLPXLA18bio.mat phy zoo detri z_rho 

figure
subplot(1,3,1),...
plot(squeeze(phy(1,:)),z_rho,'r.-',squeeze(zoo(1,:)),z_rho,'b.-',squeeze(detri(1,:)),z_rho,'magenta.-')
hold on
load('S020matNB2I.mat','S0');
zarray = .25:.5:299.75;
grid on
xlim([0 .7])
set(gca,'XTick',[0 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7])
ylim([-160 -60])
set(gca,'YTick',linspace(-160,-60,11))
set(gca,'XTickLabel',{'0','','.1','','.2','','.3','','.4','','.5','','.6','','7'})
ylabel('Depth [m]','FontSize',13)
xlabel('[mmol N/m^3]','FontSize',13);
set(gca,'FontSize',13);
title('(A) ROMS','FontSize',13)

subplot(1,3,2),...
plot(squeeze(S0(601:1200)),-zarray,'r.-',squeeze(S0(1201:1800)),-zarray,'b.-',squeeze(S0(1801:2400)),-zarray,'magenta.-')
clear S0
grid on
xlim([0 .7])
set(gca,'XTick',[0 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7])
ylim([-160 -60])
set(gca,'YTick',linspace(-160,-60,11))
set(gca,'XTickLabel',{'0','','.1','','.2','','.3','','.4','','.5','','.6','','7'})
ylabel('Depth [m]','FontSize',13)
set(gca,'FontSize',13);
title('(B) 1D, 300 m domain','FontSize',13)
xlabel('[mmol N/m^3]','FontSize',13);


load S020NB2ini_dz1H600.mat S0
zarray = .5:599.5;
subplot(1,3,3),...
plot(squeeze(S0(601:1200)),-zarray,'r.-',squeeze(S0(1201:1800)),-zarray,'b.-',squeeze(S0(1801:2400)),-zarray,'magenta.-')
grid on
xlim([0 .7])
set(gca,'XTick',[0 .05 .1 .15 .2 .25 .3 .35 .4 .45 .5 .55 .6 .65 .7])
ylim([-160 -60])
set(gca,'YTick',linspace(-160,-60,11))
set(gca,'XTickLabel',{'0','','.1','','.2','','.3','','.4','','.5','','.6','','7'})
set(gca,'FontSize',13);
ylabel('Depth [m]','FontSize',13)
title('(C) 1D, 600 m domain','FontSize',13)
xlabel('[mmol N/m^3]','FontSize',13);

%% Figure S2, grid resolution test 1D model
clear all

load 1DNPZD_fullw_out.mat
load 1DNPZD_fullw_params.mat

figure
plot(tmat,squeeze(.5.*sum(nanmean(Ooutmat(601:2400,:,18),3),1))...
    ,'k-','linewidth',2)
hold on

load 1DNPZD_fullw_dz1.mat
load 1DNPZD_fullw_dz1_params.mat

plot(tmat,squeeze(sum(nanmean(Ooutmat(601:2400,:,18),3),1))...
    ,'r--','linewidth',2)

grid on
xlabel('[Days]')
ylabel('Biomass [mmol N/m^2]')
ylim([30 100])
xlim([0 72])
legend('H=300 m, \Delta z = .5m','H=600 m, \Delta z = .5m')
title('Depth-integrated biomass, y = 5 km','fontsize',13,'fontweight','normal');
set(gca,'fontsize',13)
%% Figure 11
clear all

load 1DNPZD_fullw_out.mat
load 1DNPZD_fullw_params.mat

figure
plot(tmat,squeeze(.5.*sum(nanmean(Ooutmat(601:2400,:,13:25),3),1))...
    ,'k-','linewidth',2)
hold on


% load 2D ROMS output
load COSLPXLA18bio.mat bio ocean_time z_w 
dzgw = repmat(diff(z_w,1,1),[1 600]);

plot(ocean_time./86400,sum(dzgw.*squeeze(nanmean(bio(100:141,:,:),1)),1)...
     ,'r-','linewidth',2)

clear bio

% load 2D ROMS output for CM simulation
load COSLPXHA18_NVMbio.mat bio 

plot(ocean_time./86400,sum(dzgw.*squeeze(nanmean(bio(100:141,:,:),1)),1)...
     ,'g--','linewidth',2)
 
 clear bio dzgw z_w
 

grid on
xlabel('[Days]')
ylabel('[mmol N/m^2]')
%ylim([5 45])
ylim([35 65])
xlim([0 72])
title('Average depth-integrated P_{int}+Z_{int}+D_{int}, y=0 to 12 km','fontweight','normal','fontsize',13)
legend('1D','2D','2D CM')
set(gca,'fontsize',13)

%% Figure 10
clear all
load COSLPXLA18bio.mat bio z_rho ocean_time z_w y_rho
ytars = -12:12;
y_r = y_rho./1000-60;
dzgw = repmat(diff(z_w,1,1),[1 600]);
hfroms = fspecial('average',[1 65]);


load 1DNPZD_fullw_out.mat
load 1DNPZD_fullw_params.mat


hfpowell = fspecial('average',[1 8]);


figure

for ij = 16:2:22
[y,yidx] = min(abs(y_r-ytars(ij)));
yidx
subplot(3,4,(ij-14)./2),...
plot(tmat,imfilter(.5.*sum(Ooutmat(601:2400,:,ij),1),hfpowell,'replicate','same','corr')...
    ,'k-','linewidth',2)
hold on;
 plot(ocean_time./86400,imfilter(sum(dzgw.*squeeze(nanmean(bio(yidx-1:yidx+1,:,:),1)),1),hfroms,'replicate','same','corr')...
     ,'r--','linewidth',2)
grid on
xlabel('[Days]')
ylabel('8-d filtered biomass [mmol N/m^2]')
ylim([30 100])
xlim([0 72])
title(strcat(num2str(ytars(ij)),'km'))
end


load 1DNPZD_meanw_out.mat

for ij = 16:2:22
subplot(3,4,ij./2-7),...
plot(tmat,imfilter(.5.*sum(Ooutmat(601:2400,:,ij),1),hfpowell,'replicate','same','corr')...
    ,'g--','linewidth',2)
hold on;
end


load 1DNPZD_pertw_out.mat

for ij = 16:2:22
subplot(3,4,ij./2-7),...
hfpowell = fspecial('average',[1 8]);
plot(tmat,imfilter(.5.*sum(Ooutmat(601:2400,:,ij),1),hfpowell,'replicate','same','corr')...
    ,'g-','linewidth',2)
hold on;
end

% load ROMS vertical velocities
load w_filt.mat

for ij = 16:2:22
    subplot(3,4,4+(ij-14)./2),...
hold on
     plot(oceantime,zarray(:,ij)...
     ,'r-','linewidth',2)
hold on
    hfromsw = fspecial('average',[194 1]);
 plot(oceantime,imfilter(zarray(:,ij),hfromsw,'replicate','same','corr')...
     ,'r--','linewidth',2)
ylabel('Depth [m]')
grid on
xlabel('[Days]')
xlim([0 72.71])
ylim([-10 25])
    subplot(3,4,8+(ij-14)./2),...
hold on
     plot(oceantime,warray(:,ij)...
     ,'r-','linewidth',2)
hold on
    hfroms = fspecial('average',[194 1]);
 plot(oceantime,imfilter(warray(:,ij),hfroms,'replicate','same','corr')...
     ,'r--','linewidth',2)
ylabel('Vertical Velocity [m/s]')
grid on
xlabel('[Days]')
xlim([0 72.71])
ylim([-15 15])
end

