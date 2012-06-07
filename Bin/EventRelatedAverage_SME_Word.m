%% NOTE: Have to run this twice, setting the index by hand


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITIONS - Paths and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SbjList = {'lx2_s1_b8_mcmst_sk','lx2_s2_b8_mcmst','lx2_s3_b8_mcmst','lx2_s4_b8_mcmst', 'lx2_s5_b8_sss', 'lx02_s06_b8_mc','lx2_s07_b8_mc','lx02_s08_b8_mc', 'lx2_s09_b8_mc', 'lx02_s10_b8_mc', 'lx02_s11_b8_mc', 'lx02_s12_b8_mc'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load each sbj
for sbj = 1:16,

clear dat data d index1 index2 index_lf index_lr index_eu index_es ;

cfg = [];
%cfg.filename = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/', SbjList{sbj}, '_wrd.mat');
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
    
% Language: This is zeroing out the indexed part, so it ends up 1=ES, 2=EU
index(:,lngidx) = logical(index(:,lngidx).*0); % *2,*4 = ES; *1,*3 = EU

% Append the data
cfg = [];
dat = appenddata(cfg, d{index});

cfg = [];
cfg.channel     = {'MEG*'};
if sbj==14,
    cfg.bpfilter      = 'yes';
    cfg.bpfreq        = [0.5 10];
    cfg.bpfiltdir     = 'twopass';
elseif sbj==2 || sbj==9
    cfg.bpfilter      = 'yes';
    cfg.bpfreq        = [0.5 40];
    cfg.bpfiltdir     = 'twopass';
end;

dat = ft_preprocessing(cfg, dat);


% NOTE: For sbj 7, we excluded trials at this point using ft_rejectvisual
%       But for rejection-series-2, it was not necessary


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

% Exclude rejected trials from indexes
index1 = n1';

% Collapse and average
cfg = [];
cfg.channel     = {'MEG*'};
cfg.trials      = index1;
t{sbj,lngidx} = ft_timelockanalysis(cfg, dat); 

if lngidx==2,
    ga{sbj} = ft_timelockgrandaverage(cfg, t{sbj,1:2}); 
end;

% Apply baseline
cfg = [];
cfg.baseline = [-0.2 0.0];
t{sbj,lngidx} = ft_timelockbaseline(cfg, t{sbj,lngidx} ); 

if lngidx==2,
    ga{sbj} = ft_timelockbaseline(cfg, ga{sbj}); 
end;

end;

% % Read in this template average file, using fiff_read_evoked_all
% if sbj<=9
%     fname = strcat('Ave/lx0', int2str(sbj),'_wrd.fif'); 
% else
%     fname = strcat('Ave/lx', int2str(sbj),'_wrd.fif'); 
% end;
% 
% [data] = fiff_read_evoked_all(fname);    
% 
% % Copy the fieldtrip averages to the template 
% % 
% if sbj <= 5
%     data.evoked(1).epochs(1:306,:) = t{sbj,1}.avg;
%     data.evoked(2).epochs(1:306,:) = t{sbj,2}.avg;
%     data.evoked(3).epochs(1:306,:) = ga{sbj}.avg;
% else
%     data.evoked(1).epochs(1:306,:) = t{sbj,1}.avg(1:306,1:299);
%     data.evoked(2).epochs(1:306,:) = t{sbj,2}.avg(1:306,1:299);
%     data.evoked(3).epochs(1:306,:) = ga{sbj}.avg(1:306,1:299);
% end;
% 
% 
% data.evoked(1).comment = 'ES';
% data.evoked(2).comment = 'EU';
% data.evoked(3).comment = 'Grand Average';
% 
% % Save using fiff_write_evoked
% if sbj<=9
%     fname = strcat('Ave/lx0', int2str(sbj),'_word.fif');
% else
%     fname = strcat('Ave/lx', int2str(sbj),'_word.fif');
%    
% end;
% 
% data.info.filename = fname;
% fiff_write_evoked(fname, data);


end;



% Calculate trial numbers
for sbj=1:16,
	for cnd=1:2,
		try, TrialN(sbj,cnd) = t{sbj,cnd}.dof(1); end;
	end
end




% Make grand average
cfg = [];
cfg.keepindividual = 'no';
cfg.latency = [-0.2 0.8];
t1 = ft_timelockgrandaverage(cfg, t{[1:16],1});
t2 = ft_timelockgrandaverage(cfg, t{[1:16],2});






%% Summary averages for later-remembered and later-forgotten
d
% MAG: Scale as fT and setup a difference of averages
LR = t2; LF = t1;
LR.avg = LR.avg.*1e15;
LF.avg = LF.avg.*1e15;
LR.individual = LR.individual.*1e15;
LF.individual = LF.individual.*1e15;
diff = LR;
diff.avg = LR.avg - LF.avg;
diff.individual = LR.individual-LF.individual;
diff2 = LF;
diff2.avg = LF.avg - LR.avg;
diff2.individual = LF.individual-LR.individual;


% PLANAR: Scale as fT and setup a difference of averages
LR = t2; LF = t1;
LR.avg = LR.avg.*1e13;
LF.avg = LF.avg.*1e13;
LR.var = LR.var.*1e13;
LF.var = LF.var.*1e13;
LR.individual = LR.individual.*1e13;
LF.individual = LF.individual.*1e13;
diff = LR;
diff.avg = LR.avg - LF.avg;
diff.individual = LR.individual-LF.individual;
diff2 = LF;
diff2.avg = LF.avg - LR.avg;
diff2.individual = LF.individual-LR.individual;


% Both MAG and PLANAR: Scale as fT
LR = t1; LF = t2;
LR.avg(3:3:306,:) = t1.avg(3:3:306,:).*1e15; % MAG
LF.avg(3:3:306,:) = t2.avg(3:3:306,:).*1e15;
% LR.individual(:,3:3:306,:) = t1.individual(:,3:3:306,:).*1e15;
% LF.individual(:,3:3:306,:) = t2.individual(:,3:3:306,:).*1e15;
LR.avg(1:3:306,:) = t1.avg(1:3:306,:).*1e13; % PLAN1
LF.avg(1:3:306,:) = t2.avg(1:3:306,:).*1e13;
% LR.individual(:,1:3:306,:) = t1.individual(:,1:3:306,:).*1e13;
% LF.individual(:,1:3:306,:) = t2.individual(:,1:3:306,:).*1e13;
LR.avg(2:3:306,:) = t1.avg(2:3:306,:).*1e13; % PLAN2
LF.avg(2:3:306,:) = t2.avg(2:3:306,:).*1e13;
% LR.individual(:,2:3:306,:) = t1.individual(:,2:3:306,:).*1e13;
% LF.individual(:,2:3:306,:) = t2.individual(:,2:3:306,:).*1e13;
diff = LR;
diff.avg = LR.avg - LF.avg;
% diff.individual = LR.individual-LF.individual;
diff2 = LR;
diff2.avg = LF.avg - LR.avg;
% diff2.individual = LF.individual-LR.individual;




% Butterfly, average
subplot(2,2,1)
plot(LF.time, LF.avg)
title('Euskaraz')
axis tight
ylim([-200 200])
axis square
subplot(2,2,2)
plot(LR.time, LR.avg)
title('Gaztelaniaz')
axis tight
ylim([-200 200])
axis square

% Butterfly, variance
subplot(2,2,3)
plot(log(abs(LF.avg)), log(LF.var), '.')
title('EU')
axis tight
axis square
ylabel('Log \sigma')
xlabel('Log ABS(\mu)')
subplot(2, 2,4)
plot(log(abs(LR.avg)), log(LR.var), '.')
title('GA')
axis tight
axis square




% Plot subset
figure;
cfg = [];
% cfg.channel     = {'MEG1611' 'MEG1621' 'MEG1641' 'MEG1631'};
% cfg.channel     = {'MEG*1'};
cfg.channel     = {'MEG0242' 'MEG0233' 'MEG1513' '1612' '1623' '1812' '1823'};
% cfg.channel     = {'MEG1512' 'MEG0242' 'MEG1613' 'MEG1622' '1813' '1643' '1632' '1843'};
% cfg.channel     = LR.label(find(stat.negclusterslabelmat==1));
%cfg.channel = {'MEG0223'};
cfg.graphcolor  = [0 0.75 0; 0.75 0 0];
ft_singleplotER(cfg, LR, LF); title(' ');
xlabel('s')
ylabel('fT/cm')

[legh objh outh outm] = legend('GA','EU');
set(legh, 'FontSize', 14);
set(legh, 'EdgeColor', [1 1 1]);
set(legh, 'DataAspectRatio', [.1 0.2 2])
set(legh, 'Location', 'NorthWest');
axis square; axis tight;







% CI Plot (for PLANAR)
chidx = [6    10    12    14    16   110   118   120   134];
% chidx = find(stat.negclusterslabelmat==1);
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
legend({'LR','LF'},'Location','SouthEast', 'FontSize', FntSz);
legend boxoff;
 





% Plot individual subject
figure;
for sbj=1:16,
    subplot(4,4,sbj);
    
    cfg = [];
    cfg.channel     = {'MEG1611' 'MEG1621' 'MEG1641' 'MEG1631'};
    % cfg.channel     = {'MEG1512' 'MEG0242' 'MEG1613' 'MEG1622' '1813' '1643' '1632' '1843'};
    % cfg.channel     = {'MEG0242' 'MEG0233' 'MEG1513' '1612' '1623' '1812' '1823'};
    % cfg.channel   = LR.label(find(stat.negclusterslabelmat==1));
    

    cfg.graphcolor  = [0.85 0.15 0; 0 0.5 0];
    cfg.zparam        = 'avg';
    
    ft_singleplotER(cfg, t{sbj,1}, t{sbj,2});
    
    axis square; axis tight;
    %ylim([-0.25e-11 1.25e-11])
    title(strcat('s',int2str(sbj)));

    if sbj==1,
        xlabel('s');
        ylabel('T/cm');
        [legh objh outh outm] = legend('GA','EU');
        set(legh, 'FontSize', 14);
        set(legh, 'EdgeColor', [1 1 1]);
        set(legh, 'DataAspectRatio', [.1 0.2 2])
        set(legh, 'Location', 'NorthWest');
    end;
end;




% Layout plot
h = figure;
cfg = [];
%cfg.layout      = 'Config/NM306all.CBS.lay';
%cfg.layout      = 'Config/NM306cmb.lay';
cfg.layout      = 'Config/NM306mag.CBS.lay';
cfg.layout      = 'Config/NM306planar-lon.lay';
cfg.xparam      = 'time';
cfg.yparam      = 'avg';
% cfg.ylim        = [-1e-11 1e-11];
cfg.xlim        = [-0.2 0.8];
% cfg.linewidth   = 1;
% cfg.graphcolor  = [0.85 0.15 0; 0 0.5 0];
cfg.showlabels  = 'yes';

ft_multiplotER(cfg, t{9,1}, t{9,2});

% s14 has problematic lo-frq noise
figure;
ft_multiplotER(cfg, ga{9:16});



% Layout plot
h = figure;
cfg = [];
cfg.layout      = 'Config/NM306all.CBS.lay';
cfg.layout      = 'Config/NM306planar-lon.lay';
% cfg.layout      = 'Config/NM306mag.CBS.lay';
cfg.layout      = 'Config/NM306planar.CBS.lay';
cfg.xparam      = 'time';
cfg.yparam      = 'avg';
cfg.ylim        = 'maxmin';
cfg.linewidth   = 1;
cfg.graphcolor  = [0.85 0.15 0; 0 0.5 0];
cfg.showlabels   = 'yes';
cfg.axes        = 'no';

EU = LF; 
GA = LR;
ft_multiplotER(cfg, GA, EU);





% Topo (quick, see below for detailed version)
cfg = [];
cfg.layout      = 'Config/NM306planar-lon.lay';
% cfg.layout      = 'Config/NM306mag.CBS.lay';
cfg.zparam      = 'avg';
cfg.zlim        = [-200 200];
cfg.zlim        = 'maxmin';
cfg.xparam      = 'time';
% cfg.xlim        = [0.11 0.17];
cfg.xlim        = [0.2 0.8];
cfg.style      = 'both';
% figure;
subplot(2,2,3); ft_topoplotER(cfg, LF); title('EU'); colorbar; axis equal; axis square;
subplot(2,2,4); ft_topoplotER(cfg, LR); title('ES'); colorbar; axis equal; axis square;










%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistics 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.channel             = {'MEG*2'};
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
% cfg.layout              ='Config/NM306mag.CBS.lay';
cfg.neighbourdist       = 0.25;

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

[stat] = ft_timelockstatistics(cfg, t2, t1);









cfg = [];
cfg.alpha  = 0.05;
cfg.zparam = 'stat';
% cfg.zparam = 'ref';
cfg.zlim   = 'maxmin';
cfg.layout = 'Config/NM306planar-lon.lay';
% cfg.layout ='Config/NM306mag.CBS.lay';
ft_clusterplot(cfg, stat);



















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simple means 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


diff2.time(140:225)
cfg.timeindex = 140:225;
laterDiff = squeeze(mean(mean(diff2.individual(:, cfg.channel, cfg.timeindex),3),2));





