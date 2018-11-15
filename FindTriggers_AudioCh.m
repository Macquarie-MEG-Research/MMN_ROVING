% path            = '/Volumes/White_WD_Pa/Child_MMN/MMN_good/2632/ReTHM';
% file            = '2632_SA_ME125_2017_08_08_B1_denoise_rethm.con';
% Threshold       = 0.25;
% Sampling_Rate   = 1000;
% First_Channel   = 146; % First trigger channel
% Last_Channel    = 158; % Last trigger channel
% audio_Channel   = 135; % -1 for none
% audio_Threshold = 0.25; % Photo detector or Audio sensor threshold
% cfg.dataset = file;
% hdr.Fs=Sampling_Rate;
% cfg.trialdef.prestim=0.1;
% cfg.trialdef.poststim=0.5;


% cfg                     = [];
% cfg.dataset             = filename;
% cfg.path                = pathname;
% cfg.trialdef.prestim    = 1;         % pre-stimulus interval
% cfg.trialdef.poststim   = 1;        % post-stimulus interval
% cfg.trialdef.eventtype  = 'trigger';
% cfg.trialdef.eventvalue = 1:7;%Trigger numbers
% cfg.trialfun            = 'FindTriggers_AudioCh';
%cfg.Fs
% cfg                     = ft_definetrial(cfg);
%cfg.fixed_offset; for no audio

function [trl, event] = FindTriggers_AudioCh(cfg)

Audio_Threshold = 0.25; %safe to hardcode this as gets binarised
Threshold       = 2.0;

%function [trl, event] = FindTriggers_AudioCh(path,file,Fs,First_Channel,Last_Channel,audio_Channel,prestim,poststim)

ffile   = fullfile(cfg.path, cfg.dataset);
fprintf(1,'FindTriggers: Processing %s\n', ffile);
data    = getYkgwData(ffile);
samples = data(cfg.First_Channel:cfg.Last_Channel,:)';

[numsamples,numchannels] = size(samples);

if (cfg.Audio_Channel ~= 0)
    detector = data(cfg.Audio_Channel,:)';
    
    
    [filt_detector,ylower] = envelope(detector,50,'rms'); %calc envelope of detector
    
    [a,b]   = peakfinder(filt_detector,.05);
    yinterp = interp1(a,b,1:length(filt_detector),'nearest');
    
    bfilt_detector = filt_detector>0.5*yinterp'; %binarise detector at 50% max threshold
    timestamp      = find(diff([-1;bfilt_detector;-1]) ~= 0); % where does V change
    runlength      = diff(timestamp) ;
    timestamp0     = find(diff([-1;bfilt_detector;-1]) > 0); % where does V change
    runlength0     = runlength(1+(bfilt_detector(1)==1):2:end);
    
    for i =1:length(runlength0)-1
        if runlength0(i)<200
            bfilt_detector(timestamp0(i):timestamp0(i+1)) = 1; %gets rid of short intervals
        end
    end
    
    detector = bfilt_detector;
end

events = [];
er     = 1; % er holds the current row in the events matrix

for s = 2:numsamples
    for c = 1:numchannels
        if samples(s,c) > Threshold && samples(s-1,c) <= Threshold
            % A trigger occurred
            time = s;
            
            % Correct for photo detector if a audio channel was given,
            % but don't correct the audio channel itself!
            if (cfg.Audio_Channel ~= 0 && c ~= (cfg.Audio_Channel - cfg.First_Channel + 1))
                % If the photo detector is high when the trigger
                % occurred wait until it goes low before starting
                % the wait for it to go high
                while time <= numel(detector) && detector(time) > Audio_Threshold
                    time = time + 1;
                end
                
                % Wait for the photo detector to go above the
                % threshold
                while time <= numel(detector) && detector(time) <= Audio_Threshold
                    time = time + 1;
                end
                
                if time > numel(detector)
                    fprintf('FindTriggers: Warning - A trigger occurs at Tsec %6.3f but the audio channel does not exceed the threshold (%g) at or after this point\n', (s/cfg.Fs), Audio_Threshold);
                elseif (time - s) > 75
                    fprintf('FindTriggers: Warning - Photo detector event for trigger at Tsec %6.3f is more than 75ms later at %6.3f.\n', (s/cfg.Fs), (time/cfg.Fs));
                end
            end
            
            
            events(er,:) = [(time/cfg.Fs) c (time - s)];
            er           = er + 1;
        end
    end
end

events = sortrows(events);

if cfg.Audio_Channel == 0 %no audio
    events(:,1)=events(:,1)+cfg.fixed_offset/cfg.Fs; %add offset time to correct for delay
end

delays = events(:,3);

fprintf('FindTriggers: Output -%6.0f audio detector events of%6.0f were deemed to be outliers.\nReplaced by the median delay%6.0fms'...
    , (length(find(delays>median(delays)+3*mad(delays)))+length(find(delays<median(delays)-3*mad(delays)))), length(events), median(delays));

%triggers3 = load(['oddball_short_',cfg.dataset(1:4),'_run1.txt']);%load textfile with trigger coding


event6    = repmat(struct('type','trigger','value',1,'sample',1,'time',1,'duration','duration'), 1, length(events));

rr=diff(events(:,2));
c=1;
for i=1:length(rr)
    if rr(i)==0
        c=c+1;
        triggers_n(i)=c;
        
    else
        c=1;
        triggers_n(i)=c;
    end
end

triggers_n=[1 triggers_n]; %first is lost by diff
% 
% 
% if  ~(length(triggers3) == length(events))
%     fprintf('FindTriggers: Warning - %6.3f too many events found.\nWill try to fix it up', (length(events)-length(triggers3)))
%     if (length(events)-length(triggers3)) > 50
%         error('Big discrepancy between triggers recorded and triggers in .txt file. Stopping here')
%     elseif (length(events)-length(triggers3)) < 0
%         error('Triggers are missing from recording. Stopping here')
%     end
% end
% 
% while ~(length(triggers3) == length(events))
%     erridx = find(triggers3(:,3)-events(1:length(triggers3),2),1);
%     events = [events(1:erridx-1,:);events(erridx+1:end,:)];
% end

%events(find(events(:,3)>median(delays)+3*mad(delays) | events(:,3)<median(delays)-3*mad(delays)),3) = median(delays);

%find outlying delays - likely caused by poor detection on audio channel
%and correct by using the media,

correction_long  = events(find(events(:,3)>median(delays)+3*mad(delays)))-median(delays);
correction_short = events(find(events(:,3)<median(delays)-3*mad(delays)))-median(delays);
idx_long         = find(events(:,3)>median(delays)+3*mad(delays));
idx_short        = find(events(:,3)<median(delays)-3*mad(delays));

for t= 1:length(correction_long)
    events(idx_long(t),1) = events(idx_long(t),1)-correction_long(t)/cfg.Fs;
end

for t= 1:length(correction_short)
    events(idx_short(t),1) = events(idx_short(t),1)-correction_short(t)/cfg.Fs;
end

for k = 1:length(triggers_n)
    event6(k).type   = 'trigger';
    event6(k).value  = triggers_n(k);
    event6(k).sample = events(k,1) * cfg.Fs;
    event6(k).time   = events(k,1);
end

event = event6;

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * cfg.Fs);
posttrig = round(cfg.trialdef.poststim * cfg.Fs);
sample   = round([event(find(strcmp('duration', {event.duration}))).sample]');


trl = [];
for j = 1:length(triggers_n)
    %trg1     = value(n);
    trlbegin = sample(j) + pretrig;
    trlend   = sample(j) + posttrig;
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
end

end
