function [trl, event] = trig_fun_ME125(cfg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function for Fieldtrip trial selection was created for ARC Discovery Project (DP170103148) - 
% MEG rovling MMN study (MEG_ID: 125). It can be easily adapted for use for other MEG data acquaried
% at KIT-Macquarie Brain Reserach Liboratory. 
% Sub-function that creates an .evt file was written by Steven Saunders. 
%
% Wei He & Paul Sowman
% 2018-08-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = cfg.dataset;
path = cfg.path;
First_Channel = 146; % Frist trigger channel
Last_Channel = 158; % Last trigger channel
PD_Channel = 135; % -1 for none 
PD_Threshold = 0.25; % Photo detector or Audio sensor threshould

hdr   = ft_read_header(cfg.dataset,'dataformat','yokogawa_con');
triggers     = ft_read_data(cfg.dataset,'dataformat','yokogawa_con','chanindx',First_Channel:Last_Channel);
for i=1:size(triggers,1)
            trig_height(i)=max(triggers(i,:));
end
trig_thresh=0.25*max(trig_height);


% Threshold = 2; 
Threshold = 0.5*max(trig_height); 
Sampling_Rate = hdr.Fs; 

FindTriggers(path, file, First_Channel,Last_Channel,Threshold,PD_Channel,PD_Threshold,Sampling_Rate)

event_ori = ft_read_event(cfg.dataset,'dataformat','yokogawa_con','trigindx',First_Channel:Last_Channel,'threshold',trig_thresh,'detectflank','up');

% eventfile=[cfg.dataset(1:end-4),'.evt'];
event5 = importdata([cfg.dataset(1:end-4),'.evt']);
event6 = repmat(struct('type','trigger','value',1,'sample',1,'time',1,'duration','duration'), 1, length(event_ori));

triggers3 = load(['oddball_short_',cfg.dataset(1:4),'_run1.txt']);%load textfile with trigger coding

for k = 1:length(event_ori)
    event6(k).type = event_ori(k).type;
    event6(k).value = triggers3(k,4);
    event6(k).sample = event5.data(k,1) * hdr.Fs;
    event6(k).time = event5.data(k,1);
end

event = event6;

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig = round(cfg.trialdef.poststim * hdr.Fs);
sample = round([event(find(strcmp('duration', {event.duration}))).sample]');


trl = [];
for j = 1:length(event_ori)
    %trg1     = value(n);
    trlbegin = sample(j) + pretrig;
    trlend   = sample(j) + posttrig;
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subfunctions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function FindTriggers(path, file, First_Channel,Last_Channel,Threshold,PD_Channel,PD_Threshold, Sampling_Rate)
% Export the triggers from Yokogawa ".con" files into tab-delimited .evt files

% Change log:
%
% 2015-05-27: If the photo detector is above the threshold when a trigger
% occurrs then wait until it drops below the threshold and then change the
% time of the trigger to the time when the photo detector goes back up
% above the threshold
% 2014-05-05: Make the sampling rate user-adjustable

        ffile = fullfile(path, file);
        fprintf(1,'FindTriggers: Processing %s\n', ffile);
        data = getYkgwData(ffile);
        samples=data(First_Channel:Last_Channel,:)'; 
        
        [numsamples,numchannels] = size(samples);
       
        if (PD_Channel ~= 0)
            detector=data(PD_Channel,:)';
        end
        
        events = [];
        er = 1; % er holds the current row in the events matrix
        
        for s = 2:numsamples,
            for c = 1:numchannels,  
                if samples(s,c) > Threshold && samples(s-1,c) <= Threshold
                    % A trigger occurred
                    time = s;
                    
                    % Correct for photo detector if a PD channel was given,
                    % but don't correct the PD channel itself!
                    if (PD_Channel ~= 0 && c ~= (PD_Channel - First_Channel + 1))
                        % If the photo detector is high when the trigger
                        % occurred wait until it goes low before starting
                        % the wait for it to go high
                        while time <= numel(detector) && detector(time) > PD_Threshold
                            time = time + 1;
                        end
                        
                        % Wait for the photo detector to go above the
                        % threshold
                        while time <= numel(detector) && detector(time) <= PD_Threshold
                            time = time + 1;
                        end
                            
                        if time > numel(detector),
                            fprintf('FindTriggers: Warning - A trigger occurs at Tsec %6.3f but the PD channel does not exceed the threshold (%g) at or after this point\n', (s/Sampling_Rate), PD_Threshold);
                        elseif (time - s) > 55,
                            fprintf('FindTriggers: Warning - Photo detector event for trigger at Tsec %6.3f is more than 55ms later at %6.3f.\n', (s/Sampling_Rate), (time/Sampling_Rate));
                        end
                    end
                    
                    events(er,:) = [(time/Sampling_Rate) c 1];
                    er = er + 1;
                end
            end
        end        
        
        events = sortrows(events);
        
        % Construct output filename
        [~,name,~] = fileparts(file);
        outfname = fullfile(path,[name '.evt']);
        
        if exist(outfname) ~= 0,
            button = questdlg(['"' outfname '" already exists. Do you want to replace it?'], 'FindTriggers', 'Cancel', 'Replace', 'Cancel');
            
            if strcmp(button, 'Cancel')
                fprintf('FindTriggers: Aborted - file %s already exists\n', outfname);
                return;
            else
                delete(outfname);
            end
        end
        
        fprintf(1,'FindTriggers: Creating %s\n', outfname);
        
        fid=fopen(outfname, 'w');
        fprintf(fid, 'Tsec \t TriNo \t Code \n');
        fclose(fid);
        
        dlmwrite(outfname, events, 'delimiter', '\t', 'precision', '%6.3f', '-append');
        
    end

end

