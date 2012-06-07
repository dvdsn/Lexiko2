function [ dat ] = ftrv( dat )
% FTRV Call to ft_rejectvisual with cfg filled in
%   See ft_rejectvisual for explanation of cfg

cfg = [];
cfg.channel     = {'MEG*2' 'MEG*3'};
cfg.method      = 'summary';
cfg.keepchannel = 'yes';
cfg.metric      = 'max';

dat = ft_rejectvisual(cfg, dat);

end

