% This function computes the phase stability index, or phase stationarity
% index based on Mazaheri et al. 200X PNAS
% It uses a windowed approach to assess the similarity of phase between two time points
% t1 and t1+ncycles; ncycles is the period length of the frequency of interest;
% if an oscillation is perfectly stationary then the phase at t1 and t1+ncycles will be identical
% leading to a resultant vector length = 1; 
% therefore the higher the phase stability index is, the more stationary
% the oscillation is; 

function [pstb,pstbz]=Phase_Stab2(cfg,data)

f=cfg.freq;
width=cfg.width;
sr=cfg.sr;
ncycles = cfg.ncyc;
time=data.time;
zfoi=nearest(f,40);
zfoi=[1:zfoi-10 zfoi+10:numel(f)];

cnt=0;
for freq = f
    cnt=cnt+1;
    freq
    % decompose the signal
    cfg= []; cfg.output = 'fourier'; cfg.method = 'wavelet';
    cfg.keeptrials = 'yes'; cfg.width = width; cfg.foi =freq;
    cfg.toi = 'all';
    % get the time freqeuency decomposition
    tfr = ft_freqanalysis(cfg, data);
    
    phase_length = round(sr/freq)*ncycles;
    
    %normalize the fspctrm and cut useless information
    
    fspctrm = single(squeeze(tfr.fourierspctrm(:,:,:, nearest(tfr.time, tfr.time(1)+1): nearest(tfr.time, tfr.time(end)-1))));
    clear tfr;
    fspctrm = fspctrm./abs(fspctrm);
    
    time_axis = nan(size(fspctrm));
    nan_ind = find(squeeze(~isnan(fspctrm(1,1,:))))';
    fspctrm = reshape(fspctrm((~isnan(fspctrm))),size(fspctrm,1), []);
    
    tmp = nan(size(fspctrm));
    for tim = 1 : 1: size(fspctrm, 1)-phase_length
        
        phase_irreg = abs(mean(fspctrm([tim tim+phase_length]),1));
        
        tmp(tim+round(phase_length/2))=phase_irreg;
    end
    
    pstb(1,cnt)=nanmean(tmp,1);
end

% Apply z-transformation where mean and std are calculated outside a
% frequency of interest (i.e. 40 Hz)
zmn=mean(pstb(1,zfoi),2);
zstd=std(pstb(1,zfoi),0,2);
pstbz=(pstb-zmn)./zstd;