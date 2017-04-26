% generate time-filtered vertical velocity for 1D NPZD model
% from raw ROMS float output
% Whitt, Taylor and Levy (2017) Synoptic to planetary scale wind variab
% ility enhances phytoplankton biomass at ocean fronts
% J. Geophys. Res. Oceans

clear all
clear global
filtlen = 48; % 48 hour boxcar filter window
hf = fspecial('average',[1 filtlen]);
%ztar = 104.5;


load LF_flts_raw.mat  % load raw ROMS float trajectories
ot=repmat(oceantime,[1 length(yflt(:,1))])';

warray = zeros(1800,25);
zarray = zeros(1800,25);
yarray = -12:12;
yparray = zeros(1800,25);

ytars = -12:12;
for ytar = ytars
    zfltnow = zflt;
    yfltnow = yflt;
    mask2 = logical(abs(yflt(:,5)-ytar)>.5);
    zfltnow(mask2,:) = nan;
    yfltnow(mask2,:) = nan;
    sum(isnan(zfltnow(:,5)) == 0)
    zfltfilt = squeeze(imfilter(nanmean(zfltnow,1)-nanmean(zflt(:,5),1),hf,'replicate','same','corr'));
    yfltfilt = squeeze(imfilter(nanmean(yfltnow,1)-nanmean(yflt(:,5),1),hf,'replicate','same','corr'));
    warray(2:end,ytar+13) = diff(zfltfilt,1,2)./diff(oceantime,1,1)';
    zarray(1:end,ytar+13) = zfltfilt;
    yparray(1:end,ytar+13) = yfltfilt;
    clear zfltfilt zfltnow yfltnow mask2
    pause(.2)
end

save w_filt.mat oceantime warray yarray yparray zarray