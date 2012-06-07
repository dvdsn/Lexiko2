%% NOTE: Have to run this twice, setting the index by hand


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITIONS - Paths and variables (if necessary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load each sbj
for sbj = 1:16,

clear dat data d index1 index2 index_lf index_lr index_eu index_es ;

cfg = [];
% cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/', SbjList{sbj}, '_wrd.mat');
cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/lx02_s', int2str(sbj),'_pwr_ret_art-rejection-series-2.mat');

load(cfg.filename); 
clear dat index n1;
    
% Restore trialdef and index
% Restore trialdef
for blk=1:9,
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
cfg.channel     = {'MEG*'};
% if sbj==14,
%     cfg.bpfilter      = 'yes';
%     cfg.bpfreq        = [0.5 10];
%     cfg.bpfiltdir     = 'twopass';
% elseif sbj==2 || sbj==9
%     cfg.bpfilter      = 'yes';
%     cfg.bpfreq        = [0.5 40];
%     cfg.bpfiltdir     = 'twopass';
% end;
dat = ft_preprocessing(cfg, dat);

% Simple threshold artifact rejection (NOTE: Either grad XOR mag, not both)
p = single(cell2mat(dat.trial));                          % data
p_md = max([median(p(1:3:306,:)) median(p(2:3:306,:))]);  % md
p_sd = max([std(p(1:3:306,:)) std(p(2:3:306,:))]);        % sd
art_threshold = p_md + 6.*(p_sd);                         % 6
clear p p_md p_sd;

% Make an index for rejected trials
for k=1:size(dat.trial,2),
    if sum(max(abs([dat.trial{k}(1:3:306,:)' dat.trial{k}(2:3:306,:)']))>art_threshold) > 1
        n1(k) = false;
    else
        n1(k) = true;
    end;
end;

% Pick trials
index_lf = dat.trialinfo(:,2) == 201;   % Response to later-forgotten
index_lr = dat.trialinfo(:,2) == 202;   % Response to later-remembered

index_es = logical(rem(dat.trialinfo(:,3),2));        % *7 (MEG rsp to ES-wrd)
index_eu = logical(abs(rem(dat.trialinfo(:,3),2)-1)); % *8 (MEG rsp to EU-wrd)

index1 = index_lf & index_eu;
index2 = index_lr & index_eu;
index3 = index_lf | index_lr | index_es;

% Exclude rejected trials from indexes
index1 = index1 & n1';
index2 = index2 & n1';
index3 = index3 & n1';


% Collapse and average
cfg = [];
cfg.channel     = {'MEG*2' 'MEG*3'};
% cfg.channel     = {'MEG*1'};
cfg.method     = 'wltconvol';
cfg.output     = 'pow';
cfg.taper      = 'hanning';
cfg.foi        = [3:47];
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.5;
% cfg.tapsmofrq  = 
cfg.toi        = -2.0:0.0125:2.0;

cfg.trials     = index1;
f{sbj,1} = ft_freqanalysis(cfg, dat);

cfg.trials     = index2;
f{sbj,2} = ft_freqanalysis(cfg, dat);

end;







% Calculate trial numbers doesnt work yet
for sbj=1:16,
	for cnd=1:2,
		try, TrialN(sbj,cnd) = f{sbj,cnd}.dof(1); end;
	end
end





% Calculate baseline
cfg = [];
cfg.baseline        = [-0.9 -0.1];
cfg.baselinetype    = 'absolute';
cfg.param           = 'powspctrm';
for k=1:size(f,1),
        F{k,1} = ft_freqbaseline(cfg, f{k,1});
        F{k,2} = ft_freqbaseline(cfg, f{k,2});
end;



% Calculate grand average over subjects
cfg = [];
cfg.keepindividual = 'no';
% cfg.foilim         = [3:47];
% cfg.toilim         = [-0.25 1.25];
cfg.channel        = {'MEG*1'};
cfg.channel          = {'MEG*2' 'MEG*3'};

f1 = ft_freqgrandaverage(cfg, F{[1:10 12:16],1});
f2 = ft_freqgrandaverage(cfg, F{[1:10 12:16],2});

diff = f2;
diff.powspctrm = f1.powspctrm - f2.powspctrm;
diff.grad = f0.grad;
diffp = ft_combineplanar([], diff);

f0 = ft_freqgrandaverage(cfg, F{[1:16],1:2});
f0.grad = f{1,1}.grad;
f0p = ft_combineplanar([], f0);



%% Plotting


% grand average
figure
subplot(1,2,1)
imagesc(f1.time, f1.freq, squeeze(mean(f0p.powspctrm,1)));
axis xy
subplot(1,2,2)
imagesc(f2.time, f2.freq, squeeze(mean(f2.powspctrm,1)));
axis xy




% Simple plot
cfg = [];
% cfg.baseline     = [-0.9 -0.1]; 
% cfg.baselinetype = 'absolute'; 
% cfg.zlim         = [-1.5e-25 1.5e-25];	
cfg.zlim         = [-5e-22 5e-22];	
%cfg.zlim         = [-0.5 0.5];
cfg.xlim         = [-0.1 1.5];
cfg.ylim         = [3 40];
cfg.showlabels   = 'no';	
cfg.layout       = 'Config/NM306all.CBS.lay';
% cfg.layout       = 'Config/NM306mag.CBS.lay';
cfg.layout      = 'Config/NM306planar-lon.lay';
cfg.layout      = 'Config/NM306cmb.lay';
cfg.box          = 'yes';
figure 
ft_multiplotTFR(cfg, diffp);





% Single average
cfg = [];	
% cfg.baseline     = [-0.9 -0.1]; 
% cfg.baselinetype = 'relchange'; 
cfg.zlim         = [-2e-22 2e-22];	
% cfg.zlim         = [-0.25 0.25];	
cfg.channel      = {'MEG0242+0243' 'MEG0232+0233' 'MEG0212+0213' }; % Theta
% cfg.channel      = {'MEG2043' 'MEG1922' 'MEG2112' }; % Alpha, beta
cfg.xlim         = [-1.0 1.0];
figure 
ft_singleplotTFR(cfg, diffp);



% Topoplot
figure;
cfg = [];
cfg.baseline     = [-0.9 -0.1]; 
cfg.baselinetype = 'absolute'; 
cfg.xlim         = [1.0 1.5];   
cfg.zlim         = [-2.0e-22 2.0e-22];
% cfg.zlim         = [-0.5 0.5];
cfg.ylim         = [8 14];
cfg.layout       = 'Config/NM306mag.CBS.lay';
cfg.layout      = 'Config/NM306planar-lon.lay';
ft_topoplotTFR(cfg, f0); title('');



figure
subplot(1,2,1)
cfg.layout      = 'Config/NM306planar-lat.lay';
ft_topoplotTFR(cfg, f0); title('lat');
subplot(1,2,2)
cfg.layout      = 'Config/NM306planar-lon.lay';
ft_topoplotTFR(cfg, f0); title('lon');

figure
cfg.layout      = 'Config/NM306cmb.lay';
ft_topoplotTFR(cfg, f0p); title('cmb');






% Combined plot

subplot(2,2,1)
cfg = [];	
cfg.baseline     = [-0.9 -0.1]; 
cfg.baselinetype = 'absolute'; 
cfg.zlim         = [-2e-22 2e-22];	
cfg.channel      = {'MEG0242+0242' 'MEG0232+0233' 'MEG0212+0213' }; % Theta
cfg.xlim         = [-0.1 1.5];
cfg.masknans      = 'yes';
ft_singleplotTFR(cfg, f0p); axis square;
title('Left Frontal');
xlabel('s');
ylabel('Hz');

subplot(2,2,2)
cfg = [];
cfg.baseline     = [-0.9 -0.1]; 
cfg.baselinetype = 'absolute'; 
cfg.xlim         = [0.1 0.5];   
cfg.zlim         = [-2.0e-22 2.0e-22];
cfg.ylim         = [3 7];
cfg.layout       = 'Config/NM306cmb.lay';
cfg.comment      = 'xlim';
ft_topoplotTFR(cfg, f0p); title('\theta (3-7 Hz)');


subplot(2,2,3)
cfg = [];	
cfg.baseline     = [-0.9 -0.1]; 
cfg.baselinetype = 'absolute'; 
cfg.zlim         = [-2e-22 2e-22];	
cfg.channel      = {'MEG2042+2043' 'MEG1922+1923' 'MEG2112+2113' }; % Alpha, beta
cfg.xlim         = [-0.1 1.5];
cfg.masknans      = 'yes';
ft_singleplotTFR(cfg, f0p);  axis square;
title('Posterior');
xlabel('s');
ylabel('Hz');

subplot(2,2,4)
cfg = [];
cfg.baseline     = [-0.9 -0.1]; 
cfg.baselinetype = 'absolute'; 
cfg.xlim         = [0.5 1.5];   
cfg.zlim         = [-2.0e-22 2.0e-22];
cfg.ylim         = [8 12];
cfg.layout      = 'Config/NM306cmb.lay';
cfg.comment      = 'xlim';
ft_topoplotTFR(cfg, f0p); title('\alpha (8-12 Hz)');










% Plot individual subject
figure;
for sbj=1:16,
    subplot(4,4,sbj);
    
    cfg = [];
    cfg.baseline     = [-0.9 -0.1]; 
    cfg.baselinetype = 'absolute'; 
    cfg.zlim         = [-1e-21 1e-21];      
    cfg.xlim         = [-0.1 1.5];
    cfg.ylim         = [3 20];
%     cfg.channel      = {'MEG0242' 'MEG0233' 'MEG0212' }; % Theta
    cfg.channel      = {'MEG2042+2043' 'MEG1922+1923' 'MEG2112+2113' }
    
    diff = f{sbj,1};
    diff.powspctrm = f{sbj,2}.powspctrm-f{sbj,1}.powspctrm;
    diff.grad = f0.grad;
    diffp = ft_combineplanar([], diff);

    ft_singleplotTFR(cfg, diffp);
    
    axis square; axis tight;
    %ylim([-0.25e-11 1.25e-11])
    title(strcat('s',int2str(sbj)));

end;








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate baseline
cfg = [];
cfg.baseline        = [-0.9 -0.1];
cfg.baselinetype    = 'relative';
cfg.param           = 'powspctrm';
for k=1:size(f,1),
    
        F{k,1} = ft_freqbaseline(cfg, f{k,1});
        F{k,2} = ft_freqbaseline(cfg, f{k,2});
        
        F{k,1}.grad = f0.grad;
        F{k,2}.grad = f0.grad;

        F{k,1} = ft_combineplanar([], F{k,1});       
        F{k,2} = ft_combineplanar([], F{k,2});
        
end;


% Calculate grand average over subjects
cfg = [];
cfg.keepindividual = 'yes';
% cfg.foilim         = [3 40];
% cfg.toilim         = [-0.0 1.5];
% cfg.channel        = {'MEG*3'};
% cfg.channel          = {'MEG*2' 'MEG*3'};

f1 = ft_freqgrandaverage(cfg, F{[1:10 12:16],1});
f2 = ft_freqgrandaverage(cfg, F{[1:10 12:16],2});

f1.grad = F{1,1}.grad;
f2.grad = F{1,1}.grad;


% Cluster analysis
cfg = [];
cfg.channel             = {'MEG*3'};
cfg.latency             = [0.0 1.5];
cfg.frequency           = [3 47];
cfg.avgovertime         = 'no';
cfg.avgoverfreq         = 'no';
cfg.avgoverchan         = 'no';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 4;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.numrandomization    = 1500;
cfg.layout              = 'Config/NM306cmb.lay';
cfg.neighbourdist       = 0.25;

subj = 15;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

[stat] = ft_freqstatistics(cfg, f2, f1);

