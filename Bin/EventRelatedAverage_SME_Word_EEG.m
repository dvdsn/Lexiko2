
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITIONS - Paths and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SbjList = {'lx2_s1_b8_mcmst_sk','lx2_s2_b8_mcmst','lx2_s3_b8_mcmst','lx2_s4_b8_mcmst', 'lx2_s5_b8_sss', 'lx02_s06_b8_mc','lx2_s07_b8_mc','lx02_s08_b8_mc', 'lx2_s09_b8_mc', 'lx02_s10_b8_mc', 'lx02_s11_b8_mc'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load each sbj
for sbj = 1:16,

clear dat data index d index1 index2 index_eu index_es n1;

cfg = [];
% cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/', SbjList{sbj}, '_wrd.mat');
cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/lx02_s', int2str(sbj),'_wrd_art-rejection-series-2.mat');

load(cfg.filename); 
    
for lngidx = 1:2,

clear dat index n1;
    
% Restore trialdef and index
for blk=1:8,
    for cnd=1:2,
		try,
           		index(blk,cnd) = isstruct(d{blk,cnd});
        		d{blk,cnd}.trialdef = d{blk,cnd}.cfg.trialdef;
                d{blk,cnd}.trialdef = d{blk,cnd}.cfg.previous.trialdef;

		end;
	end;
end;
    
% Language: lngidx 1=ES, 2=EU
index(:,lngidx) = logical(index(:,lngidx).*0); % *2,*4 = ES; *1,*3 = EU

% Append the data
cfg = [];
dat = appenddata(cfg, d{index});

cfg = [];
cfg.channel     = {'EEG*'};
cfg.lpfreq      = 15;
cfg.lpfilter    = 'yes';
cfg.reref       = 'yes';
cfg.refchannel  = 'EEG010';
dat = ft_preprocessing(cfg, dat);

% Collapse and average
cfg = [];
cfg.channel     = {'EEG00*' 'EEG010'};
t{sbj,lngidx} = ft_timelockanalysis(cfg, dat); 

end;

cfg = [];
cfg.baseline = [-0.2 0.0];
t{sbj,1} = ft_timelockbaseline(cfg, t{sbj,1} ); 
t{sbj,2} = ft_timelockbaseline(cfg, t{sbj,2} );

end;



% Calculate trial numbers
for sbj=1:16,
	for cnd=1:2,
		TrialN(sbj,cnd) = t{sbj,cnd}.dof(1);
	end
end


% Make grand average
% NOTE: Sbj 5 excluded
cfg = [];
cfg.keepindividual = 'no';
cfg.latency = [-0.2 0.8];
t1 = ft_timelockgrandaverage(cfg, t{[1:4 6:16],1});
t2 = ft_timelockgrandaverage(cfg, t{[1:4 6:16],2});


% Simple plot over subjects
ChannelNames = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
figure;
for k=1:9,
    
    subplot(3,3,k)
    plot(t1.time, t1.avg(k,:), 'g');
    hold on;
    plot(t2.time, t2.avg(k,:), 'r');    
     axis([-0.2 0.8 -5e-6 5e-6]);
    axis square;
    title(ChannelNames{k});
    
    if k==1
        legend('ES','EU');
        legend boxoff;
    end;
    
end;


% Simple plot per subject
cfg = [];
cfg.channel = 4;

figure;
for k=1:16,
    
    subplot(4,4,k)
    plot(t{k,1}.time, t{k,1}.avg(cfg.channel,:), 'g');
    hold on;
    plot(t{k,2}.time, t{k,2}.avg(cfg.channel,:), 'r');    
    axis tight;
    axis square;
    title(k);
    
end;


