function [trl, event] = trig_fun_125_mmn(cfg)

trigger_nums = 146:155;
triggers     = ft_read_data(cfg.dataset,'dataformat','yokogawa_con','chanindx',trigger_nums);
%triggers(33,:) = 0; %this channel is always high and not used for triggers
%triggers(7,:)  = 0; %this channel is the audio and will create loads of triggers
hdr          = ft_read_header(cfg.dataset,'dataformat','yokogawa_con');

for j=1:size(triggers,1)
    trig_height(j) = max(triggers(j,:)); %Y = prctile(x,42)
end

trig_thresh = 0.25*max(trig_height);
triggers    = triggers>trig_thresh; %Binarise trigger channels

event = [];
list  = trigger_nums;

for k=1:size(triggers,1)
    channel   = list(k);
    trig      = triggers(k,:);
    pad       = trig(1);
    trigshift = 0;
    begsample = 1;
    
    
    for l=find(diff([pad trig(:)'])>0)
        event(end+1).type   = num2str(channel);
        event(end  ).sample = l + begsample - 1;      % assign the sample at which the trigger has gone down
        event(end  ).value  = trig(l+trigshift);      % assign the trigger value just _after_ going up
    end
    
end

triggers2 = load( 'oddball_short_2712_run1.txt');%load textfile with trigger coding
values    = triggers2(1:length(triggers2),3);

%% Fill up the event struct
event2 = nestedSortStruct(event,'sample');
event3 = {event2.type};
event3 = cell2mat(cellfun(@str2num,event3,'un',0))-126;
event3 = arrayfun(@num2str, event3, 'unif', 0);
event4 = repmat(struct('type','trigger','value',1,'sample',1,'time',1,'duration',''), 1, length(event3));

times2 = [event2.sample];

for m=1:length(event3)
    event4(m).time   = times2(m)/hdr.Fs; %%****FIXME****depends on Fs
    event4(m).sample = times2(m);
    event4(m).value  = triggers2(m,4);
end

event = event4;

% search for "trigger" events
value  = {event(find(strcmp('trigger', {event.type}))).value}';
sample = round([event(find(strcmp('trigger', {event.type}))).sample]');

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig = round(cfg.trialdef.poststim * hdr.Fs);

trl = [];
for n = 1:length(value)
    %trg1     = value(n);
    trlbegin = sample(n) + pretrig;
    trlend   = sample(n) + posttrig;
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
end

end

