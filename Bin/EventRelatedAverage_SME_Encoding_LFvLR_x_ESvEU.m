%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SME Encoding 1: Analysis of encoding trials, pooled over languages
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITIONS - Paths and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SbjList = {'lx02_s01_b8_mc_lowcor','lx2_s2_b8_mcmst','lx2_s3_b8_mcmst','lx2_s4_b8_mcmst', 'lx2_s5_b8_sss', 'lx02_s06_b8_mc'};



LngList = {'EU','ES'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load each sbj
for sbj = 7:7,

clear d data index index1 index2;

cfg = [];
%cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/', SbjList{sbj}, '_enc.mat');
cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/lx02_s', int2str(sbj), '_enc_art-rejection-series-2.mat');


load(cfg.filename); 


for lngidx = 1:1,

clear dat index n1;
    

% Restore trialdef
for blk=1:9,
    for enclng=1:4,
    for cnd=1:2,
	for rsp=1:2,
		try
           		index(blk,enclng,cnd,rsp) = isstruct(d{blk,enclng,cnd,rsp});
        		d{blk,enclng,cnd,rsp}.trialdef = d{blk,enclng,cnd,rsp}.cfg.previous.trialdef;

		end;
	end;
    end;
    end;
end;
    

% Language
index(:,lngidx+0,:,:) = logical(index(:,lngidx+0,:,:).*0); % *2,*4 = ES; *1,*3 = EU
index(:,lngidx+2,:,:) = logical(index(:,lngidx+2,:,:).*0); % *2,*4 = ES; *1,*3 = EU

% Append the data
cfg = [];
data = ft_appenddata(cfg, d{index});

% Baseline correct single trials (and apply other sbj-specif. options)
cfg = [];
cfg.channel        = {'MEG*'};
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.2 0.0];

if sbj==14,
    cfg.bpfilter      = 'yes';
    cfg.bpfreq        = [0.5 10];
    cfg.bpfiltdir     = 'twopass';
elseif sbj==2 || sbj==9
    cfg.bpfilter      = 'yes';
    cfg.bpfreq        = [0.5 40];
    cfg.bpfiltdir     = 'twopass';    
end;

data = ft_preprocessing(cfg, data);


% Simple threshold artifact rejection (NOTE: Either grad XOR mag, not both)
p = single(cell2mat(data.trial));                          % data
p_md = max([median(p(1:3:306,:)) median(p(2:3:306,:))]);  % md
p_sd = max([std(p(1:3:306,:)) std(p(2:3:306,:))]);        % sd
art_threshold = p_md + 6.*(p_sd);                         % 6
clear p p_md p_sd;

% Make an index for rejected trials
for k=1:size(data.trial,2),
    if sum(max(abs([data.trial{k}(1:3:306,:)' data.trial{k}(2:3:306,:)']))>art_threshold) > 1
        n1(k) = false;
    else
        n1(k) = true;
    end;
end;


% Index for later-rem v later-forgotten
index1 = data.trialinfo(:,3) == 201; % Response to later-forgotten
index2 = data.trialinfo(:,3) == 202; % Response to later-remembered

% Exclude rejected trials from indexes
index1 = index1 & n1';
index2 = index2 & n1';

% Collapse and average
cfg = [];
% cfg.channel     = {'MEG*2' 'MEG*3'};
cfg.channel     = {'MEG*'};
cfg.trials      = index1;
t{sbj,1} = ft_timelockanalysis(cfg, data); 
cfg.trials      = index2;
t{sbj,2} = ft_timelockanalysis(cfg, data);

ga{sbj} = ft_timelockgrandaverage(cfg, t{sbj,1:2}); 


% Baseline (if not earlier)
cfg = [];
cfg.baseline = [-0.2 0.0];
t{sbj,1} = ft_timelockbaseline(cfg, t{sbj,1} ); 
t{sbj,2} = ft_timelockbaseline(cfg, t{sbj,2} );

ga{sbj} = ft_timelockbaseline(cfg, ga{sbj}); 


% % Read in this template average file, using fiff_read_evoked_all
% fname = strcat('Ave/', 'template.fif'); 
% [data] = fiff_read_evoked_all(fname);
% 
% % Copy the fieldtrip averages to the template 
% % 
% data.evoked(1).epochs(1:306,:) = t{sbj,1}.avg;
% data.evoked(2).epochs(1:306,:) = t{sbj,2}.avg;
% data.evoked(3).epochs(1:306,:) = ga{sbj}.avg; 
% 
% data.evoked(1).comment = 'LF';
% data.evoked(2).comment = 'LR';
% data.evoked(3).comment = 'Grand Average';
% 
% % Save using fiff_write_evoked
% fname = strcat('Ave/', SbjList{sbj}(1:7),'enc_LFvLR_',LngList{lngidx},'.fif');
% data.info.filename = fname;
% fiff_write_evoked(fname, data);
% 


end;

end;



%% Average

% Notes:
% 0. Take care to track artifact rejection, esp. EOG and EMG.
% 1. Sbj 11 has very few LR, so we have to exclude the data.
% 2. Sbj 14 has (very) hi-amp >10 Hz noise artifact, so data was selectively
%    filtered with the lo-pass corner at 10 Hz, and above 0.5 Hz.



% Make grand average
cfg = [];
cfg.keepindividual = 'no';
cfg.latency = [-0.2 0.8];
cfg.subjectindex = [1:10 12:16];

t1 = ft_timelockgrandaverage(cfg, t{cfg.subjectindex,1});
t2 = ft_timelockgrandaverage(cfg, t{cfg.subjectindex,2});




% Get trial numbers
for sbj=1:16,
   for cnd=1:2,
    TrialCount(sbj,cnd) = t{sbj,cnd}.dof(1);    
   end
end






%% Summary averages for later-remembered and later-forgotten

% MAG: Scale as fT and setup a difference of averages
LR = t2; LF = t1;
LR.avg = LR.avg.*1e15;
LF.avg = LF.avg.*1e15;
diff = LR;
diff.avg = LR.avg - LF.avg;
diff2 = LR;
diff2.avg = LF.avg - LR.avg;
diff2.individual = LF.individual-LR.individual;
LR.individual = LR.individual.*1e15;
LF.individual = LF.individual.*1e15;
diff.individual = LR.individual-LF.individual;


% PLANAR: Scale as fT and setup a difference of averages
LR = t2; LF = t1;
LR.avg = LR.avg.*1e13;
LF.avg = LF.avg.*1e13;
diff = LR;
diff.avg = LR.avg - LF.avg;
% diff.individual = LR.individual-LF.individual;
diff2 = LR;
diff2.avg = LF.avg - LR.avg;
% diff2.individual = LF.individual-LR.individual;
% LR.individual = LR.individual.*1e13;
% LF.individual = LF.individual.*1e13;


% Both MAG and PLANAR: Scale as fT
LR = t2; LF = t1;
LR.avg(3:3:306,:) = t2.avg(3:3:306,:).*1e15; % MAG
LF.avg(3:3:306,:) = t1.avg(3:3:306,:).*1e15;
% LR.individual(:,3:3:306,:) = t1.individual(:,3:3:306,:).*1e15;
% LF.individual(:,3:3:306,:) = t2.individual(:,3:3:306,:).*1e15;
LR.avg(1:3:306,:) = t2.avg(1:3:306,:).*1e13; % PLAN1
LF.avg(1:3:306,:) = t1.avg(1:3:306,:).*1e13;
% LR.individual(:,1:3:306,:) = t1.individual(:,1:3:306,:).*1e13;
% LF.individual(:,1:3:306,:) = t2.individual(:,1:3:306,:).*1e13;
LR.avg(2:3:306,:) = t2.avg(2:3:306,:).*1e13; % PLAN2
LF.avg(2:3:306,:) = t1.avg(2:3:306,:).*1e13;
% LR.individual(:,2:3:306,:) = t1.individual(:,2:3:306,:).*1e13;
% LF.individual(:,2:3:306,:) = t2.individual(:,2:3:306,:).*1e13;
diff = LR;
diff.avg = LR.avg - LF.avg;
% diff.individual = LR.individual-LF.individual;
diff2 = LR;
diff2.avg = LF.avg - LR.avg;
% diff2.individual = LF.individual-LR.individual;





%% Plotting

% Butterfly, average
subplot(2,2,1)
plot(LF.time, LF.avg)
title('Later-Forgotten')
axis tight
ylim([-250 250]./1.0)
axis square
subplot(2,2,2)
plot(LR.time, LR.avg)
title('Later-Recalled')
axis tight
ylim([-250 250]./1.0)
axis square

% Butterfly, variance
subplot(2,2,3)
plot(LF.avg, log(LF.var), '.')
title('LF')
axis tight
axis square
subplot(2,2,4)
plot(LR.avg, log(LR.var), '.')
title('LR')
axis tight
axis square






% Plot subset
figure;
cfg = [];
cfg.channel     = {'MEG1611' 'MEG1621' 'MEG1641' 'MEG1631'};
cfg.channel     = {'MEG0242' 'MEG0233' 'MEG1513' '1612' '1623' '1812' '1823'};
% cfg.channel     = {'MEG1512' 'MEG0242' 'MEG1613' 'MEG1622' '1813' '1643' '1632' '1843'};
cfg.graphcolor    = [0 0.75 0; 0.75 0 0];
ft_singleplotER(cfg, LR, LF); title(' ');
xlabel('s')
ylabel('fT/cm')

[legh objh outh outm] = legend('LR','LF');
set(legh, 'FontSize', 14);
set(legh, 'EdgeColor', [1 1 1]);
set(legh, 'DataAspectRatio', [.1 0.2 2])
set(legh, 'Location', 'NorthWest');
axis square; axis tight;






% CI Plot
chidx = [14 16 110];
N = size(LR.individual,1);
FntSz = 20;

LR.se = squeeze(mean(std(LR.individual(:,chidx,:),0,1)./sqrt(N)))';
LF.se = squeeze(mean(std(LF.individual(:,chidx,:),0,1)./sqrt(N)))';

LR.mn = mean(LR.avg(chidx,:));
LF.mn = mean(LF.avg(chidx,:));

addpath Bin;
figure;
hold on;
ciplot(LR.mn+LR.se, LR.mn-LR.se, LR.time, [0.00 0.75 0.00]);
alpha(0.25);
ciplot(LF.mn+LF.se, LF.mn-LF.se, LF.time, [0.75 0.00 0.00]);
alpha(0.5);

plot(LF.time, LF.mn, 'Color', [0.75 0.0  0.0], 'LineWidth', 3)
plot(LR.time, LR.mn, 'Color', [0.0  0.75 0.0], 'LineWidth', 3)

hold off;
axis tight;
xlabel('s', 'FontSize', FntSz); ylabel('fT/cm', 'FontSize', FntSz);
legend({'LR','LF'},'Location','NorthWest', 'FontSize', FntSz);
legend boxoff;






% Plot individual subject
figure;
for sbj=1:16,
    subplot(4,4,sbj);
    
    cfg = [];
    cfg.channel     = {'MEG1611' 'MEG1621' 'MEG1641' 'MEG1631'};
    cfg.channel     = {'MEG1512' 'MEG0242' 'MEG1613' 'MEG1622' '1813' '1643' '1632' '1843'};
    cfg.channel     = {'MEG0242' 'MEG0233' 'MEG1513' '1612' '1623' '1812' '1823'};
    cfg.channel     = {'MEG1413'};
    cfg.graphcolor      = [0 0.75 0; 0.75 0 0];
    cfg.zparam        = 'avg';
    
    ft_singleplotER(cfg, t{sbj,2}, t{sbj,1});
    
    axis square; axis tight;
    title(strcat('s',int2str(sbj)));

    if sbj==1,
        xlabel('s');
        ylabel('T/cm');
        [legh objh outh outm] = legend('LR','LF');
        set(legh, 'FontSize', 14);
        set(legh, 'EdgeColor', [1 1 1]);
        set(legh, 'DataAspectRatio', [.1 0.2 2])
        set(legh, 'Location', 'NorthWest');
    end;
end;






% Layout plot
h = figure;
cfg = [];
cfg.layout      = 'Config/NM306all.CBS.lay';
cfg.layout      = 'Config/NM306planar-lon.lay';
% cfg.layout      = 'Config/NM306mag.CBS.lay';
% cfg.layout      = 'Config/NM306cmb.lay';
cfg.xparam      = 'time';
cfg.yparam      = 'avg';
cfg.ylim        = 'maxmin';
cfg.linewidth   = 1;
cfg.graphcolor  = [0.85 0.15 0; 0 0.5 0];
cfg.showlabels  = 'yes';

ft_multiplotER(cfg, LF, LR);







% Layout plot, indy subject
h = figure;
cfg = [];
cfg.layout      = 'Config/NM306all.CBS.lay';
% cfg.layout      = 'Config/NM306cmb.lay';
cfg.xparam      = 'time';
cfg.yparam      = 'avg';
cfg.ylim        = 'maxmin';
cfg.linewidth   = 1;
cfg.graphcolor  = [0.85 0.15 0; 0 0.5 0];
cfg.showlabels   = 'yes';

ft_multiplotER(cfg, t{7,1}, t{7,2}); 






% Topo
cfg = [];
cfg.layout      = 'Config/NM306planar-lon.lay';
% cfg.layout      = 'Config/NM306mag.CBS.lay';
cfg.zparam      = 'avg';
cfg.zlim        = [-50 50];
cfg.zlim        = 'maxmin';
cfg.xparam      = 'time';
cfg.xlim        = [0.4 0.8];
cfg.style      = 'both';

figure;
subplot(2,2,1); ft_topoplotER(cfg, LF); title('LF'); colorbar; axis equal; axis square; 
subplot(2,2,2); ft_topoplotER(cfg, LR); title('LR'); colorbar; axis equal; axis square;
cfg.zlim = [-50 50];
subplot(2,2,3); ft_topoplotER(cfg, diff2); title('LF-LR'); colorbar; axis equal; axis square; 










%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.channel             = {'MEG*'};
cfg.latency             = [0.0 0.8];
cfg.avgovertime         = 'no';
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
cfg.layout              ='Config/NM306planar-lon.lay';
cfg.layout              ='Config/NM306mag.CBS.lay';
cfg.neighbourdist       = 0.15;

subj = 22;
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

[stat] = ft_timelockstatistics(cfg, LR, LF);







cfg = [];
cfg.alpha  = 0.01;
cfg.zparam = 'stat';
% cfg.zparam = 'ref';
cfg.zlim   = [-4 4];
cfg.layout = 'Config/NM306planar-lon.lay';
% cfg.layout ='Config/NM306cmb.lay';
cfg.layout              ='Config/NM306mag.CBS.lay';
ft_clusterplot(cfg, stat);









% Alt plot of sig clusters (NOT WORKING YET)
figure;  
timestep        = 0.05;      % (in seconds)
sampling_rate   = LR.fsample;
sample_count    = length(stat.time);

j = [0:timestep:1];
m = [1:timestep*sampling_rate:sample_count];

pos_cluster_pvals   = [stat.posclusters(:).prob];
pos_signif_clust    = find(pos_cluster_pvals < stat.cfg.alpha);
pos                 = ismember(stat.posclusterslabelmat, pos_signif_clust);

neg_cluster_pvals   = [stat.negclusters(:).prob];
neg_signif_clust    = find(neg_cluster_pvals < stat.cfg.alpha);
neg                 = ismember(stat.negclusterslabelmat, neg_signif_clust);

for k = 1:11;
     subplot(4,3,k);   
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   
     cfg.zlim = [-1.0e-13 1.0e-13];   
     cfg.zlim = 'maxmin';
     neg_int = all(neg(:, m(k):m(k+1)), 2);
     cfg.highlight = 'on';
     cfg.highlightchannel = find(neg_int);       
     cfg.comment = 'xlim';   
     cfg.commentpos = 'title';   
     cfg.layout = 'Config/NM306planar-lon.lay';
     ft_topoplotER(cfg, diff);
end  





%% Plot for Conference(s)


% Topo
cfg = [];
cfg.layout      = 'Config/NM306planar-lon.lay';
%cfg.layout      = 'Config/NM306mag.CBS.lay';
cfg.zparam      = 'avg';
cfg.zlim        = [-400 400];
cfg.xparam      = 'time';
cfg.xlim        = [0.4 0.8];
cfg.style       = 'both';
cfg.marker      = 'off';
cfg.comment     = 'no';

figure;
subplot(2,2,1); ft_topoplotER(cfg, LR); title('\bf A: Later Recalled (LR)'); colorbar; axis equal; axis square;
subplot(2,2,2); ft_topoplotER(cfg, LF); title('\bf B: Later Forgotten (LF)'); colorbar; axis equal; axis square; 

cfg.highlight          = 'on';
cfg.highlightchannel   = find(stat.posclusterslabelmat==1);
cfg.highlightsymbol    = '*';
cfg.highlightcolor     = [1 1 1];
cfg.highlightsize      = 8;
cfg.highlightfontsize  = 10;

cfg.zlim = [-150 150];
subplot(2,2,3); ft_topoplotER(cfg, diff2); title('\bf C: LR-LF'); colorbar; axis equal; axis square; 

% CI Plot
subplot(2,2,4); 
chidx = cfg.highlightchannel';
N = size(LR.individual,1);
FntSz = 12;

LR.se = squeeze(mean(std(LR.individual(:,chidx,:),0,1)./sqrt(N)))';
LF.se = squeeze(mean(std(LF.individual(:,chidx,:),0,1)./sqrt(N)))';

LR.mn = mean(LR.avg(chidx,:));
LF.mn = mean(LF.avg(chidx,:));

hold on;
ciplot(LR.mn+LR.se, LR.mn-LR.se, LR.time, [0.00 0.75 0.00]);
alpha(0.25);
ciplot(LF.mn+LF.se, LF.mn-LF.se, LF.time, [0.75 0.00 0.00]);
alpha(0.5);
plot(LF.time, LF.mn, 'Color', [0.75 0.0  0.0], 'LineWidth', 3)
plot(LR.time, LR.mn, 'Color', [0.0  0.75 0.0], 'LineWidth', 3)
hold off;
axis tight; axis square;
xlabel('s', 'FontSize', FntSz); ylabel('fT/cm', 'FontSize', FntSz);
title('\bf D: Average ERF for Positive Cluster (*)');
legend({'LR','LF'},'Location','SouthEast', 'FontSize', FntSz);
legend boxoff;
 





