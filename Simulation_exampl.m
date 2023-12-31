addpath('U:\Castles\Simon\AssociativeMemoryFormation\fred4simon\tbx\fieldtrip-master')
ft_defaults

% let's make a signal that has noise for the first half, and a stationary
% sine wave for the second half
sr=1000;
t=0:1/sr:60;
midpnt=(numel(t)-1)/2;
noise1=(rand(1,midpnt)-0.5).*6;
noise2=(rand(1,numel(t(midpnt+1:end)))-0.5).*5;
sig=[noise1 sin(2*pi*t(midpnt+1:end)*40)+noise2];
figure;plot(t,sig);

% now we split the signal in two and organize it in a way that fieldtrip
% understands; data1 has noise, data2 has stationary sine waves. Note that
% the data has no trials;
dum                     = [];
dum.fsample             = sr;
dum.label               = {'dumChan1'};

data1=dum;
data2=dum;
data1.trial{1,1} = sig(1:midpnt);
data1.time{1,1} = t(1:midpnt);
data2.trial{1,1} = sig(midpnt+1:end-1);
data2.time{1,1} = t(midpnt+1:end-1);

% now let's calculate phase stability for the two datasets;
% set your parameters for the calculation
cfg.freq=20:1:60;
cfg.width=6;
cfg.sr=sr;

[~,pstb1]=Phase_Stab(cfg,data1);
[~,pstb2]=Phase_Stab(cfg,data2);

figure;plot(cfg.freq,[pstb1;pstb2]);

cfg.ncyc=3;
[pstb1r,~]=Phase_Stab2(cfg,data1);
[pstb2r,~]=Phase_Stab2(cfg,data2);

figure;plot(cfg.freq,[pstb1r;pstb2r]);
