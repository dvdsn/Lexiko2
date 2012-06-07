
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITIONS - Paths and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SbjList = {'lx2_s1_b8_mcmst_sk','lx2_s2_b8_mcmst','lx2_s3_b8_mcmst','lx2_s4_b8_mcmst', 'lx2_s5_b8_sss', 'lx02_s06_b8_mc','lx2_s07_b8_mc','lx02_s08_b8_mc', 'lx2_s09_b8_mc', 'lx02_s10_b8_mc', 'lx02_s11_b8_mc', 'lx02_s12_b8_mc'};

% SbjList = {'lx02_s01_b8_mc_lowcor','lx2_s2_b8_mcmst','lx2_s3_b8_mcmst','lx2_s4_b8_mcmst', 'lx2_s5_b8_sss', 'lx02_s06_b8_mc','lx2_s07_b8_mc','lx02_s08_b8_mc', 'lx2_s09_b8_mc', 'lx02_s10_b8_mc', 'lx02_s11_b8_mc', 'lx02_s12_b8_mc'};

SbjList = {'lx02_s01_b8_mc_lowcor',...
    'lx2_s2_b8_mcmst',...
    'lx2_s3_b8_mcmst',...
    'lx2_s4_b8_mcmst',...
    'lx2_s5_b8_sss',...
    'lx02_s06_b8_mc',...
    'lx2_s07_b8_mc',...
    'lx02_s08_b8_mc',...
    'lx2_s09_b8_mc',...
    'lx02_s10_b8_mc',...
    'lx02_s12_b8_mc'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load each sbj
for sbj = 1:11,

clear dat data index d index1 index2 index_eu index_es n1;

cfg = [];
cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/', SbjList{sbj}, '_ret.mat');

load(cfg.filename); 


% Restore trialdef
for blk=1:8,
    for cnd=1:2,
	for rsp=1:2,
		try,
           		index(blk,cnd,rsp) = isstruct(d{blk,cnd,rsp});
        		d{blk,cnd,rsp}.trialdef = d{blk,cnd,rsp}.cfg.trialdef;
                d{blk,cnd,rsp}.trialdef = d{blk,cnd,rsp}.cfg.previous.trialdef;

		end;
	end;
    end;
end;
    

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


% Pick trials
index_lf = dat.trialinfo(:,2) == 201;   % Response to later-forgotten
index_lr = dat.trialinfo(:,2) == 202;   % Response to later-remembered

index_es = logical(rem(dat.trialinfo(:,3),2));        % *7 (MEG rsp to ES-wrd)
index_eu = logical(abs(rem(dat.trialinfo(:,3),2)-1)); % *8 (MEG rsp to EU-wrd)

index1 = index_lr & index_es;
index2 = index_lf & index_es;
index3 = index_lf | index_lr | index_es;


% Collapse and average
cfg = [];
cfg.channel     = {'EEG00*' 'EEG010'};
cfg.trials      = index1;
t{sbj,1} = ft_timelockanalysis(cfg, dat); 
cfg.trials      = index2;
t{sbj,2} = ft_timelockanalysis(cfg, dat);

cfg = [];
cfg.baseline = [-0.2 0.0];
t{sbj,1} = ft_timelockbaseline(cfg, t{sbj,1} ); 
t{sbj,2} = ft_timelockbaseline(cfg, t{sbj,2} );


end;



% Calculate trial numbers
for sbj=1:11,
	for cnd=1:2,
		TrialN(sbj,cnd) = t{sbj,cnd}.dof(1);
	end
end


% Make grand average
cfg = [];
cfg.keepindividual = 'no';
cfg.latency = [-0.2 0.8];
t1 = ft_timelockgrandaverage(cfg, t{[1:4 7:11],1});
t2 = ft_timelockgrandaverage(cfg, t{[1:4 7:11],2});


% Simple plot over subjects
ChannelNames = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4'};
for k=1:9,
    
    subplot(3,3,k)
    plot(t1.time, t1.avg(k,:), 'g');
    hold on;
    plot(t2.time, t2.avg(k,:), 'r');    
     axis([-0.2 0.8 -1e-5 1e-5]);
    axis square;
    title(ChannelNames{k});
    
    if k==1
        legend('EU','ES');
        legend boxoff;
    end;
    
end;


% Simple plot per subject
cfg = [];
cfg.channel = 3;

for k=1:12,
    
    subplot(4,3,k)
    plot(t{k,1}.time, t{k,1}.avg(cfg.channel,:), 'g');
    hold on;
    plot(t{k,2}.time, t{k,2}.avg(cfg.channel,:), 'r');    
    axis tight;
    axis square;
    title(k);
    
end;


