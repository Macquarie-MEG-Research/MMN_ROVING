
%% Function to do Naatanens optimum 1. Based on LSD Neuroimage paper by Suresh in 2017
% Does not include direction deviant
% Has freq deviant +/- 10% Volume deviant +/- 10% (based on calibration of
% standard at ~75dB, duration deviant 1/3, and gap deviant ~7ms
%use e.g. optimum_MMN(1000,0) 1K standard with no triggers



%PFS March2018

function optimum_MMN_dir(base_frequency,portoutput)

% base_frequency = 1000;
% portoutput     = 0;
portoutput     = logical(portoutput); % This is a flag to indicate whether you want triggers or not 1 or 0

%%%%%%%%%%%%%%%%%%%%
%Create the sequence
%%%%%%%%%%%%%%%%%%%%

%%%%%
%Hardcoded levels
%%%%

halfk_vol       = [0.12 0.62];
onek_vol        = [0.2 1.0];
oneandhalfk_vol = [0.15 1.0];

if base_frequency==500
    vol = halfk_vol;
elseif base_frequency==1000
    vol = onek_vol;
else
    vol = oneandhalfk_vol;
end

deviant_multiplier = 30; %this determines total deviant numbers

deviants = [];
deviants = (repmat(2:9,1,deviant_multiplier))'; %make deviant sequence 4 deviants * n but 2 types (freq and intensity) are split into up and down
%% Need 2*4 & 5
deviants = [deviants;ones(deviant_multiplier,1)*4;ones(deviant_multiplier,1)*9]; %additional 4 and 5 which are gap and duration

% Do the reshuffling:
deviants = deviants(randperm(numel(deviants))); %initial shuffle
old_idx  = unique(find(diff(deviants)==0));%find repeats
while ~isempty(old_idx) %continue until no repeats
    new_idx                       = unique(setdiff(1:length(deviants),old_idx)); %find new spots
    new_idx                       = new_idx((randi(length(new_idx),length(old_idx),1)))';
    deviants([new_idx;old_idx],:) = deviants([old_idx;new_idx],:); %swap
    old_idx                       = unique(find(diff(deviants)==0));%find repeats
end

standards = ones(length(deviants),1);
seq       = [standards'; deviants'];
seq       = seq(:)';
seq       = [ones(1,15) seq];

% seq=[];
% seq=[ones(1,50)*7;ones(1,50)*8];
% seq       = seq(:)';

%%%%%%%%%%%%%%%%%%%%
%Now do the PTB audio stuff
%%%%%%%%%%%%%%%%%%%%

% Initialize Sounddriver
InitializePsychSound(1);

% Number of channels and Frequency of the sound
nrchannels = 2;
freq       = 48000;

% How many times to we wish to play the sound
repetitions = 1;

% Length of the beep
beepLengthSecs = 0.075;

% Length of the pause between beeps
beepPauseTime = 0.425;

% Start immediately (0 = immediately)
startCue = 0;

% Should we wait for the device to really start (1 = yes)
% INFO: See help PsychPortAudio
waitForDeviceStart = 1;

% Open Psych-Audio port, with the follow arguements
% (1) [] = default sound device
% (2) 1 = sound playback only
% (3) 1 = default level of latency
% (4) Requested frequency in samples per second
% (5) 2 = stereo putput
%pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels);

pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);

% Set the volume to half for this demo

% Make a beep which we will play back to the user
standard_freq = base_frequency;
standard      = MakeBeep(standard_freq, beepLengthSecs, freq); %STANDARD

lofreq_deviant = MakeBeep(standard_freq*0.9, beepLengthSecs, freq); %LOWER FREQ DEVIANT
hifreq_deviant = MakeBeep(standard_freq*1.1, beepLengthSecs, freq); %HIGHER FREQ DEVIANT

gap_deviant = standard; %gap deviant
dur_deviant = MakeBeep(standard_freq, beepLengthSecs/3, freq); %duration deviant

dir_deviant   = MakeBeep(standard_freq, beepLengthSecs - 0.0008, freq);
dir_deviant_f = MakeBeep(standard_freq, beepLengthSecs, freq);


b = 125;e=125; %e is duration

t      = 1:length(standard);
window = ((1+sin(pi*(t-b)/2/e))/2.*(t>b-e)-1 ) .*(t<=b+e)+ 1;
window = window.*fliplr(window);

t            = 1:length(dur_deviant);
short_window = ((1+sin(pi*(t-b)/2/e))/2.*(t>b-e)-1 ) .*(t<=b+e)+ 1;
short_window = short_window.*fliplr(short_window);

t          = 1:length(standard)/2-3.5*48; %
gap_window = ((1+sin(pi*(t-b)/2/e))/2.*(t>b-e)-1 ) .*(t<=b+e)+ 1;
gap_window = [gap_window zeros(1,7*48) 0 gap_window]; %pad with 1 zero insert 6ms of zeros 48000/1000*3
gap_window = gap_window.*fliplr(gap_window);

t          = 1:length(dir_deviant); %
dir_window = ((1+sin(pi*(t-b)/2/e))/2.*(t>b-e)-1 ) .*(t<=b+e)+ 1;
dir_window = dir_window.*fliplr(dir_window);

t           = 1:length(dir_deviant_f); %
dir_fwindow = ((1+sin(pi*(t-b)/2/e))/2.*(t>b-e)-1 ) .*(t<=b+e)+ 1;
dir_fwindow = dir_fwindow.*fliplr(dir_fwindow);

%need double of these - coded in sequence
gap_deviant = standard.*gap_window; % ****FIXME*****
%gap_deviant=[zeros(1000,1)' gap_deviant zeros(1000,1)'];

standard = standard.*window;

%frequency deviants
lofreq_deviant = lofreq_deviant.*window;
hifreq_deviant = hifreq_deviant.*window;

%need double of these - coded in sequence
gap_deviant = standard.*gap_window; % ****FIXME*****
%gap_deviant=[zeros(1000,1)' gap_deviant zeros(1000,1)'];

%need double of these - coded in sequence
dur_deviant = dur_deviant.*short_window; %
%dur_deviant=[zeros(1000,1)' dur_deviant zeros(1000,1)'];

%volume deviants need to calibrate for each frequency change later
lovol_deviant = standard;
hivol_deviant = standard; % Volume set for intensity deviant below

%volume deviants need to calibrate for each frequency change later
dir_deviant   = [zeros(39,1)' dir_deviant.*dir_window];
%right_deviant = [dir_deviant.*dir_window zeros(39,1)'];
dir_deviant_f = dir_deviant_f.*dir_fwindow;

%2 channel sound
beeps = [[standard;standard] [lofreq_deviant;lofreq_deviant] [hifreq_deviant;hifreq_deviant] [gap_deviant;gap_deviant]...
    [lovol_deviant;lovol_deviant] [hivol_deviant;hivol_deviant] [dir_deviant;dir_deviant_f]...
    [dir_deviant_f;dir_deviant] [dur_deviant;dur_deviant]];
% Fill the audio playback buffer with the audio data, doubled for stereo
% presentation

%PsychPortAudio('FillBuffer', pahandle, [beeps; beeps]); % FIXME change this dir_devaints need different chans
PsychPortAudio('FillBuffer', pahandle, beeps)
% Start audio playback
%PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);

beep_length = freq*beepLengthSecs;
start_idxs  = 1:beep_length:(8*beep_length+1);
start_idxs  = [start_idxs start_idxs(end)+beep_length/3];
end_idxs    = beep_length:beep_length:beep_length*8;
end_idxs    = [end_idxs end_idxs(end)+beep_length/3];%duration deviant is 1/3 length
volumes     = [ones(1,4)*0.5 vol(1) vol(2) 0.5 0.5 0.5]; %set the volumes INCLUDES SETTING INTENSITY DEVIANT
port_codes  = 33:42; % Use this array to assign port codes

% set up parallel port and trigger support
if portoutput
    [ioObj,address] = MQinitPP;
else
end

tic

for i =1:length(seq)
    PsychPortAudio('SetLoop',pahandle,start_idxs(seq(i)),end_idxs(seq(i)))
    PsychPortAudio('Volume', pahandle, volumes(seq(i)));
    [~, ~, ~, estStopTime] = PsychPortAudio('Stop', pahandle, 1, 1);
    
    % Compute new start time for follow-up beep, beepPauseTime after end of
    % previous one
    startCue = estStopTime + beepPauseTime;
    
    % Start audio playback
    PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
    
    code = port_codes(seq(i));
    
    if portoutput
        MQsendtrigger_PFS(ioObj,address,code, .05);
    else
    end
    
    % Wait for stop of playback
    PsychPortAudio('Stop', pahandle, 1, 1);
end

% Close the audio device
PsychPortAudio('Close', pahandle);

timeElapsed = toc;
display(['duration of block: ' num2str(timeElapsed)])

% %% Trigger function
%     function MQsendtrigger_PFS(ioObj,address,triggerline,duration)
%         
%         if ~exist('duration','var')
%             duration = .005;
%         end
%         
%         io64(ioObj,address,triggerline);
%         pause(duration)
%         io64(ioObj,address,0);
%         
%     end
% 
% %% Initialise the port Address D020 for stim 1
%     function [ioObj,address] = MQinitPP
%         
%         %create IO32 interface object
%         ioObj = io64;
%         
%         % check the port status
%         status = io64(ioObj);
%         
%         if status == 0
%             address = hex2dec('D020');
%             display(['Yay, we have a functioning parallel port at: ' num2str(address)])
%         else
%             display('Failed')
%             ioObj = [];
%         end
%     end
end