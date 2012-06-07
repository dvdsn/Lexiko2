% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESSING - Read single trial data with respect to trigger
% 
% The point of this analysis is to get data for a general word response
%                 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Notes
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% For windows
addpath G:\Lexiko\Lexiko2\Bin     % For trialfun_neuromag.m

% For linux
addpath /home/BCBL/public/Lexiko/Lexiko2/Bin

% List of subjects

SbjList = {{'lx2_s1_b1_mcmst_sk','lx2_s1_b2_mc_sk','lx2_s1_b3_mcmst_sk','lx2_s1_b4_mcmst_sk','lx2_s1_b5_mcmst_sk', 'lx2_s1_b6_mcmst_sk','lx2_s1_b7_mc_sk','lx2_s1_b8_mcmst_sk'},...
           {'lx2_s2_b1_mc','lx2_s2_b2_mcmst','lx2_s2_b3_mc','lx2_s2_b4_mcmst','lx2_s2_b5_mcmst','lx2_s2_b6_mc','lx2_s2_b7_mcmst','lx2_s2_b8_mcmst'},...
           {'lx2_s3_b1_mcmst','lx2_s3_b2_mc','lx2_s3_b3_mcmst','lx2_s3_b4_mcmst','lx2_s3_b5_mcmst','lx2_s3_b6_mcmst','lx2_s3_b7_mc','lx2_s3_b8_mcmst'},...
           {'lx2_s4_b1_mc','lx2_s4_b2_mc_lowcor','lx2_s4_b3_mcmst','lx2_s4_b4_mcmst','lx2_s4_b5_mcmst', 'lx2_s4_b6_mc_lowcor','lx2_s4_b7_mcmst','lx2_s4_b8_mcmst'},...
           {'lx2_s5_b1_sss','lx2_s5_b2_sss','lx2_s5_b3_sss','lx2_s5_b4_sss','lx2_s5_b5_sss','lx2_s5_b6_sss','lx2_s5_b7_sss','lx2_s5_b8_sss'},...
           {'lx02_s06_b1_sss','lx02_s06_b2_sss','lx02_s06_b3_mc','lx02_s06_b4_mc','lx02_s06_b5_mc','lx02_s06_b6_mc','lx02_s06_b7_mc','lx02_s06_b8_mc'},...
           {'lx02_s07_b1_mc','lx02_s07_b2_mc','lx2_s07_b3_mc','lx2_s07_b4_mcmst','lx2_s07_b5_mc','lx2_s07_b6_mc','lx2_s07_b7_mc', 'lx2_s07_b8_mc'},...
           {'lx02_s08_b1_mc','lx02_s08_b2_mc','lx02_s08_b3_mc','lx02_s08_b4_mc','lx02_s08_b5_mc','lx02_s08_b6_mc','lx02_s08_b7_mc','lx02_s08_b8_mc'},...
           {'lx2_s09_b1_mcmst','lx2_s09_b2_mc','lx2_s09_b3_mc', 'lx2_s09_b4_mc', 'lx2_s09_b5_mc', 'lx2_s09_b6_mc', 'lx2_s09_b7_mc', 'lx2_s09_b8_mc'},...
           {'lx02_s10_b1_mc','lx02_s10_b2_mc','lx02_s10_b3_mc', 'lx02_s10_b4_mc', 'lx02_s10_b5_mc', 'lx02_s10_b6_mc', 'lx02_s10_b7_mc', 'lx02_s10_b8_mc'},...
           {'lx02_s11_b1_mc','lx02_s11_b2_mc','lx02_s11_b3_mc','lx02_s11_b4_mc','lx02_s11_b5_mc','lx02_s11_b6_mc','lx02_s11_b7_mc','lx02_s11_b8_mc'},...
           {'lx02_s12_b1_mc','lx02_s12_b2_mc','lx02_s12_b3_mc','lx02_s12_b4_mc','lx02_s12_b5_mc','lx02_s12_b6_mc','lx02_s12_b7_mc','lx02_s12_b8_mc'} };


% Testing version, different artifact rejection schemes
SbjList = {{'lx02_s01_b01_mc','lx02_s01_b02_mc','lx02_s01_b03_mc','lx02_s01_b04_mc','lx02_s01_b05_mc','lx02_s01_b06_mc','lx02_s01_b07_mc','lx02_s01_b08_mc'},...
   {'lx02_s02_b01_mc','lx02_s02_b02_mc','lx02_s02_b03_mc','lx02_s02_b04_mc','lx02_s02_b05_mc','lx02_s02_b06_mc','lx02_s02_b07_mc','lx02_s02_b08_mc'},...
   {'lx02_s03_b01_mc','lx02_s03_b02_mc','lx02_s03_b03_mc','lx02_s03_b04_mc','lx02_s03_b05_mc','lx02_s03_b06_mc','lx02_s03_b07_mc','lx02_s03_b08_mc'},...
   {'lx02_s04_b01_mc','lx02_s04_b02a_mc','lx02_s04_b02b_mc','lx02_s04_b03_mc','lx02_s04_b04_mc','lx02_s04_b05_mc','lx02_s04_b06_mc','lx02_s04_b07_mc','lx02_s04_b08_mc'},...
   {'lx02_s05_b01_mc','lx02_s05_b02_mc','lx02_s05_b03_mc','lx02_s05_b04_mc','lx02_s05_b05_mc','lx02_s05_b06_mc','lx02_s05_b07_mc','lx02_s05_b08_mc'},...
   {'lx02_s06_b01_mc','lx02_s06_b02_mc','lx02_s06_b03_mc','lx02_s06_b04_mc','lx02_s06_b05_mc','lx02_s06_b06_mc','lx02_s06_b07_mc','lx02_s06_b08_mc'},...
   {'lx02_s07_b01_mc','lx02_s07_b02_mc','lx02_s07_b03_mc','lx02_s07_b04_mc','lx02_s07_b05_mc','lx02_s07_b06_mc','lx02_s07_b07_mc','lx02_s07_b08_mc'},...
   {'lx02_s08_b01_mc','lx02_s08_b02_mc','lx02_s08_b03_mc','lx02_s08_b04_mc','lx02_s08_b05_mc','lx02_s08_b06_mc','lx02_s08_b07_mc','lx02_s08_b08_mc'},...
   {'lx02_s09_b01_mc','lx02_s09_b02_mc','lx02_s09_b03_mc','lx02_s09_b04_mc','lx02_s09_b05_mc','lx02_s09_b06_mc','lx02_s09_b07_mc','lx02_s09_b08_mc'},...
   {'lx02_s10_b01_mc','lx02_s10_b02_mc','lx02_s10_b03_mc','lx02_s10_b04_mc','lx02_s10_b05_mc','lx02_s10_b06_mc','lx02_s10_b07_mc','lx02_s10_b08_mc'},...
   {'lx02_s11_b01_mc','lx02_s11_b02_mc','lx02_s11_b03_mc','lx02_s11_b04_mc','lx02_s11_b05_mc','lx02_s11_b06_mc','lx02_s11_b07_mc','lx02_s11_b08_mc'},...
   {'lx02_s12_b01_mc','lx02_s12_b02_mc','lx02_s12_b03_mc','lx02_s12_b04_mc','lx02_s12_b05_mc','lx02_s12_b06_mc','lx02_s12_b07_mc','lx02_s12_b08_mc'},...
   {'lx02_s13_b01_mc','lx02_s13_b02_mc','lx02_s13_b03_mc','lx02_s13_b04_mc','lx02_s13_b05_mc','lx02_s13_b06_mc','lx02_s13_b07_mc','lx02_s13_b08_mc'},...
   {'lx02_s14_b01_mc','lx02_s14_b02_mc','lx02_s14_b03_mc','lx02_s14_b04_mc','lx02_s14_b05_mc','lx02_s14_b06_mc','lx02_s14_b07_mc','lx02_s14_b08_mc'},...
   {'lx02_s15_b01_mc','lx02_s15_b02_mc','lx02_s15_b03_mc','lx02_s15_b04_mc','lx02_s15_b05_mc','lx02_s15_b06_1_mc','lx02_s15_b07_mc','lx02_s15_b08_mc'},...
   {'lx02_s16_b01_mc','lx02_s16_b02_mc','lx02_s16_b03_mc','lx02_s16_b04_mc','lx02_s16_b05_mc','lx02_s16_b06_mc','lx02_s16_b07_mc','lx02_s16_b08_mc'}};


% Notes on SbjList:
%
% 1. From subject s4, because we could not get enough trials, we deleted: 
%    'lx2_s4_b2_part1_mcmst' 
%    'lx2_s4_b6_mcmst_correct' 
%

       
% Key: In general, EU is odd, ES is even. Eight word-pairs per phase.
%
% *1,*3 EU in study phase
% *2,*4 ES in study phase
% *5    EU in test phase
% *6    ES in test phase
% *7    EU-prompt in test phase
% *8    ES-prompt in test phase
% 201   incorrect or no response
% 202   correct recall

ConditionList = {[11 21 31 41 51 61 71 81 13 23 33 43 53 63 73 83 15 25 35 45 55 65 75 85 ],
                 [12 22 32 42 52 62 72 82 14 24 34 44 54 64 74 84 16 26 36 46 56 66 76 86 ]};

for Sbj   = 16:16,
    
    for Block = 1:size(SbjList{Sbj},2),

        for Cnd = 1:2,
        
            
cfg = [];
cfg.channel             = {'MEG' 'EEG' 'EEG061' 'EEG062' 'EEG063' 'EEG064'};
cfg.dataset             = strcat('/home/ddavidson/ddavidson/Lexiko2/Data_MF_MC/',...
                            SbjList{Sbj}{Block}, '.fif');
                                              
cfg.dataformat          = 'neuromag_mne';
cfg.headerformat        = 'neuromag_mne';
cfg.eventformat         = 'neuromag_mne';
cfg.continuous          = 'yes';
cfg.lpfilter            = 'yes';
cfg.lpfreq              = 40;
cfg.hpfilter            = 'no';
cfg.hpfreq              = 5;
cfg.hpfilttype          = '/home/ddavidson/ddavidson/Lexiko2/Config/hipass_ret.txt';
cfg.hpfiltdir           = 'twopass';
cfg.demean              = 'no';
cfg.padding             = 3.0;

cfg.trialdef.eventvalue = ConditionList{Cnd};

% ERF:
% cfg.trialdef.prestim    = 0.2;
% cfg.trialdef.poststim   = 0.8; 

% PWR:
cfg.trialdef.prestim    = 2.0;
cfg.trialdef.poststim   = 2.0; 

cfg.trialfun           = 'trialfun_neuromag_lexiko_word';
cfg.trialdef.eventtype = 'STI101';
cfg.trialdef.trgchan   = 'STI101';

% EOG artifact
cfg.artfctdef.eog.channel         = {'EEG061' 'EEG062'};
cfg.artfctdef.eog.bpfilter        = 'yes';
cfg.artfctdef.eog.bpfilttype      = 'but';
cfg.artfctdef.eog.bpfreq          = [1 15];
cfg.artfctdef.eog.bpfiltord       = 4;
cfg.artfctdef.eog.hilbert         = 'yes';
cfg.artfctdef.eog.cutoff          = 5;
cfg.artfctdef.eog.trlpadding      = 1.5;
cfg.artfctdef.eog.fltpadding      = 1.1;
cfg.artfctdef.eog.artpadding      = 1.1;

% Jump artifact
cfg.artfctdef.jump.channel        = {'MEG*' '-MEG0143'};
cfg.artfctdef.jump.medianfilter   = 'yes';
cfg.artfctdef.jump.medianfiltord  = 7;
cfg.artfctdef.jump.absdiff        = 'yes';
cfg.artfctdef.jump.cutoff         = 20;

% Muscle (lip) artifact
cfg.artfctdef.muscle.bpfilter    = 'yes';
cfg.artfctdef.muscle.bpfreq      = [20 90];
cfg.artfctdef.muscle.bpfiltord   = 10;
cfg.artfctdef.muscle.bpfilttype  = 'but';
cfg.artfctdef.muscle.hilbert     = 'yes';
cfg.artfctdef.muscle.boxcar      = 0.2;
cfg.artfctdef.muscle.channel     = {'EEG063'};
cfg.artfctdef.muscle.cutoff      = 4;  
cfg.artfctdef.muscle.trlpadding  = 3.0;
cfg.artfctdef.muscle.fltpadding  = 3.0;
cfg.artfctdef.muscle.artpadding  = 3.0;

% Clipping artifact
cfg.artfctdef.clip.channel       = {'MEG*'};
cfg.artfctdef.clip.pretim        = 0.2;
cfg.artfctdef.clip.psttim        = 0.8;
cfg.artfctdef.clip.thresh        = 0.010;
cfg.continuous                   = 'yes';

% Rejection
cfg.artfctdef.feedback            = 'no';
cfg.artfctdef.reject              = 'complete';
cfg.artfctdef.minaccepttim        = 0.1;

cfg = ft_definetrial(cfg);

try
    
    [cfg, artifact] = ft_artifact_clip(cfg);
    %[cfg, artifact] = ft_artifact_eog(cfg);
    [cfg, artifact] = ft_artifact_jump(cfg);
    %[cfg, artifact] = ft_artifact_muscle(cfg);
    cfg = ft_rejectartifact(cfg);
    
catch
    
    cfg.trl = cfg.trl(2:end-1,:);
    [cfg, artifact] = ft_artifact_clip(cfg);
    %[cfg, artifact] = ft_artifact_eog(cfg);
    [cfg, artifact] = ft_artifact_jump(cfg);
    %[cfg, artifact] = ft_artifact_muscle(cfg);
    cfg = ft_rejectartifact(cfg);

end;

d{Block, Cnd} = ft_preprocessing(cfg); 

% % Resample
% cfg = [];
% cfg.resamplefs = 250;
% cfg.detrend    = 'no';
% d{Block, Cnd} = ft_resampledata(cfg, d{Block, Cnd});

        end;
    end;

% ERF:    
% fdat = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/lx02_s',...
%               int2str(Sbj),...
%               '_wrd_art-rejection-series-2.mat');

% Pwr          
fdat = strcat('/home/ddavidson/ddavidson/Lexiko2/Data/lx02_s',...
              int2str(Sbj),...
              '_pwr_wrd_art-rejection-series-2.mat');
          
          
save(fdat,'d','-v7.3');

fclose('all');

clear d;

end;



