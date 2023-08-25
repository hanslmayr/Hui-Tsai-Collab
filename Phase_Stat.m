% This function computes the phase stability index, or phase stationarity
% index based on Michelmann et al. 2016; PLoS Biology;
% It uses a windowed approach to assess deviation of phase from uniform
% distribution; the higher the index the less it deviates from uniform
% distribution and the more stationary the signal is;


function [pstb]=Phase_Stab(cfg,data)

freq=cfg.freq;
width=cfg.width;
sr=cfg.sr;
time=cfg.time;

for freq = 2:1:40
    freq
    % decompose the signal
    cfg= []; cfg.output = 'fourier'; cfg.method = 'wavelet';
    cfg.keeptrials = 'yes'; cfg.width = width; cfg.foi =freq;
    cfg.toi = time;
    % get the time freqeuency decomposition
    tfr = ft_freqanalysis(cfg, data);
    
    phase_length = round(sr/freq);
    
    %normalize the fspctrm and cut useless information
    
    fspctrm = single(squeeze(tfr.fourierspctrm(:,:,:, nearest(tfr.time, -1): nearest(tfr.time, 4))));
    clear tfr;
    fspctrm = fspctrm./abs(fspctrm);
    
    time_axis = nan(size(fspctrm));
    nan_ind = find(squeeze(~isnan(fspctrm(1,1,:))))';
    fspctrm = reshape(fspctrm((~isnan(fspctrm))),size(fspctrm,1),size(fspctrm,2), []);
    
    tmp = nan(size(fspctrm));
    for tim = 1 : 1: size(fspctrm, 3)-phase_length
        
        phase_irreg = abs(mean(fspctrm(:,:,tim:tim+phase_length),3));
        
        tmp(:,:,tim+round(phase_length/2))=phase_irreg;
    end
    
    time_axis(:,:,nan_ind) = tmp;
    
    hits = cat(4,...
        hits, mean(time_axis(trials_ret_mem,:,:),1)...
        );
    cr = cat(4,...
        cr, mean(time_axis(trials_ret_cr,:,:),1)...
        );
    if length(trials_ret_miss)>=15;
        miss = cat(4,...
            miss, mean(time_axis(trials_ret_miss,:,:),1)...
            );
        h_m_all = cat(4,...
            h_m_all, mean(time_axis([trials_ret_mem; trials_ret_miss],:,:),1)...
            );
    end
    hcr_all = cat(4,...
        hcr_all, mean(time_axis([trials_ret_mem; trials_ret_cr],:,:),1)...
        );
end