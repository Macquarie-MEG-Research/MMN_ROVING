%Requires psychtoolbox
%Simple function to create an abitrary duration sine at one of 3
%frequencies (500 1000 1500)Hz
%e.g. tonecal_MMN(1000,5) --> a 1k tone for 5 seconds
%PFS March2018


function tonecal_MMN(frequency,duration) 

if ~ismember(frequency,[500 1000 1500])
    error ('Illegal duration: please use 500, 1000 or 1500 only')
else
end

InitializePsychSound(1);

% Number of channels and Frequency of the sound
nrchannels = 2;
freq       = 48000;

% How many times to we wish to play the sound
repetitions = 1;

% Length of the tone
toneLengthSecs = duration; % in seconds

% Length of the pause between tones
tonePauseTime = 0.5;

% Start immediately (0 = immediately)
startCue = 0;

% Initialise PPA object
pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);

% Build the tone
tone = MakeBeep(frequency, toneLengthSecs, freq);

% Put waveform in buffer
PsychPortAudio('FillBuffer', pahandle, [tone; tone]);

% Set the calibration level this is the standard mid range
PsychPortAudio('Volume', pahandle,0.5); % dont change this volume

% Start playback
PsychPortAudio('Start', pahandle, repetitions, startCue, 1);

% Wait for stop of playback
PsychPortAudio('Stop', pahandle, 1, 1);

% Close the audio device
PsychPortAudio('Close', pahandle);

end
