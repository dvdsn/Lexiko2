
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITIONS - Paths and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SbjList = {'lx02_s01_b8_mc_lowcor',...
           'lx2_s2_b8_mcmst',...
           'lx2_s3_b8_mcmst',...
           'lx2_s4_b8_mcmst',...
           'lx2_s5_b8_sss',...
           'lx02_s06_b8_mc'};

SbjList = {'lx02_s01_b8_mc_lowcor','lx2_s2_b8_mcmst','lx2_s3_b8_mcmst','lx2_s4_b8_mcmst', 'lx2_s5_b8_sss', 'lx02_s06_b8_mc','lx2_s07_b8_mc','lx02_s08_b8_mc', 'lx2_s09_b8_mc', 'lx02_s10_b8_mc', 'lx02_s11_b8_mc'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load each sbj
for sbj = 1:11,

clear d data index index1 index2;

cfg = [];
cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/', SbjList{sbj}, '_enc.mat');

load(cfg.filename); 

 
% Restore trialdef
for blk=1:9,
    for enclng=1:4,
    for cnd=1:2,
	for rsp=1:2,
		try,
           		index(blk,enclng,cnd,rsp) = isstruct(d{blk,enclng,cnd,rsp});
        		d{blk,enclng,cnd,rsp}.trialdef = d{blk,enclng,cnd,rsp}.cfg.previous.trialdef;

		end;
	end;
    end;
    end;
end;
    

% Encoding language
index(:,1:2,:,:) = logical(index(:,1:2,:,:).*0); % *1,*3 = EU
%index(:,3:4,:,:) = logical(index(:,3:4,:,:).*0); % *2,*4 = ES

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



% *7    previous-ES-word-prompt-for-EU-rsp in test phase
% *8    previous-EU-word-prompt-for-ES-rsp in test phase
% 201   incorrect or no response
% 202   correct recall

% Pick trials

% Index for EU v ES words
%index1 = logical(abs(rem(data.trialinfo(:,1),2)));   % *7
%index2 = logical(abs(rem(data.trialinfo(:,1),2)-1)); % *8

% THis would work for retrieval, but how does it work for encoding??


% Index for later-rem v later-forgotten
index1 = dat.trialinfo(:,2) == 201; % Response to later-forgotten
index2 = dat.trialinfo(:,2) == 202; % Response to later-remembered

% Exclude rejected trials from indexes
index1 = index1;
index2 = index2;


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
for sbj=1:6,
	for cnd=1:2,
		TrialN(sbj,cnd) = t{sbj,cnd}.dof(1);
	end
end




% Make grand average
cfg = [];
cfg.keepindividual = 'no';
cfg.latency = [-0.2 0.8];
t1 = ft_timelockgrandaverage(cfg, t{[1:4 6],1});
t2 = ft_timelockgrandaverage(cfg, t{[1:4 6],2});


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
        legend('LF','LR');
        legend boxoff;
    end;
    
end;


% Simple plot per subject
cfg = [];
cfg.channel = 2;

for k=1:6,
    
    subplot(2,3,k)
    plot(t1.time, t{k,1}.avg(cfg.channel,:), 'g');
    hold on;
    plot(t2.time, t{k,2}.avg(cfg.channel,:), 'r');    
    axis tight;
    axis square;
    title(k);
    
end;


