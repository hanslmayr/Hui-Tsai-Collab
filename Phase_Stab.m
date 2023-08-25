% This function computes the phase stability index, or phase stationarity
% index based on Michelmann et al. 2016; PLoS Biology;
% It uses a windowed approach to assess deviation of phase from uniform
% distribution; the higher the index the less it deviates from uniform
% distribution and the more stationary the signal is;


function [pstb,pstbz]=Phase_Stab(cfg,data)

f=cfg.freq;
width=cfg.width;
sr=cfg.sr;
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
    
    phase_length = round(sr/freq);
    
    %normalize the fspctrm and cut useless information
    
    fspctrm = single(squeeze(tfr.fourierspctrm(:,:,:, nearest(tfr.time, tfr.time(1)+1): nearest(tfr.time, tfr.time(end)-1))));
    clear tfr;
    fspctrm = fspctrm./abs(fspctrm);
    
    time_axis = nan(size(fspctrm));
    nan_ind = find(squeeze(~isnan(fspctrm(1,1,:))))';
    fspctrm = reshape(fspctrm((~isnan(fspctrm))),size(fspctrm,1), []);
    
    tmp = nan(size(fspctrm));
    for tim = 1 : 1: size(fspctrm, 1)-phase_length
        
        phase_irreg = abs(mean(fspctrm(tim:tim+phase_length),1));
        
        tmp(tim+round(phase_length/2))=phase_irreg;
    end
    
    pstb(1,cnt)=1-nanmean(tmp,1);
end

zmn=mean(pstb(1,zfoi),2);
zstd=std(pstb(1,zfoi),0,2);
pstbz=(pstb-zmn)./zstd;