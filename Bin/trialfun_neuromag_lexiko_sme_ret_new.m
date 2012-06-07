function [trl] = trialfun_neuromag_lexiko_sme_ret_new(cfg)

% For Neuromag data, the trialdef structure should contain
%   cfg.trialdef.trgchan  = channel label, e.g. 'STI 001'
%   cfg.trialdef.prestim  = 0.300         latency in seconds
%   cfg.trialdef.poststim = 0.700         latency in seconds

% read the header information
hdr = read_header(cfg.dataset, 'headerformat', cfg.headerformat);
event = read_event(cfg.dataset, 'eventformat', cfg.eventformat);

chanindx = match_str(hdr.label, cfg.trialdef.trgchan);
nsamples = hdr.nTrials * hdr.nSamples;
checkboundary = 1;

trl = [];
retrievalblock = 0;

for k=2:length(event)-6
  if strcmp(event(k).type, 'STI101')
      
      if event(k).value == 110;
        retrievalblock = retrievalblock + 1;
      end;
      
      % For current retrieval cue and response, if conditions hold
      if sum(event(k).value == cfg.trialdef.eventvalue(:)) > 0 && ...
         sum([event((k+1):(k+6)).value] == cfg.trialdef.respvalue) > 0

      % add the probe word time index to the trl definition
      begsample = event(k-1).sample - cfg.trialdef.prestim*hdr.Fs;
      endsample = event(k-1).sample + cfg.trialdef.poststim*hdr.Fs - 1;
      offset = - cfg.trialdef.prestim*hdr.Fs;  
      
      % trl(end+1, :) = round([begsample endsample offset]);
      trl(end+1, :) = round([begsample endsample offset cfg.trialdef.respvalue ...
          cfg.trialdef.respvalue event(k).value retrievalblock]);

      end
  end
end





% Old:
% Logical index for values matching eventvalue (e.g., 21)
% trigger = dat == cfg.trialdef.eventvalue;
% trigindx = find(trigger & [0 diff(trigger)]);
% 
% 
% Read stimulus channel
% dat = read_data(cfg.dataset, ...
%     'header', hdr, ...
%     'begsample', 1, ...
%     'endsample', nsamples, ...
%     'chanindx', chanindx, ...
%     'checkboundary', checkboundary, ...
%     'dataformat', cfg.dataformat);
%
% trl(:,1) = trigindx(:) - round(cfg.trialdef.prestim*hdr.Fs);  % define begin of each trial (in samples)
% trl(:,2) = trigindx(:) + round(cfg.trialdef.poststim*hdr.Fs); % define end of each trial (in samples)
% trl(:,3) = round(-cfg.trialdef.prestim*hdr.Fs);               % define the trial offset relative to latency zero (in samples)
% trl(find(trl(:,1)<1), :) = [];                                % remove trials beginning before the start of the file
% trl(find(trl(:,2)>hdr.nTrials*hdr.nSamples), :) = [];         % remove trials after the end of the file
% fprintf('%d triggers converted into %d trials\n', length(trigindx), size(trl,1));



