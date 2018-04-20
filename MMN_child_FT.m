%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Child FT MMN is a function to preprocess MMN data collected from MQ MEG
% using roving paradigm
%
% Written by Paul Sowman Apr 2018 )
%
% INPUTS:
% - hipass        = cutoff for hipass filter

%*******************************************************
%You need 1 only .con file and the oddball run .txt file
%*******************************************************

% EXAMPLE FUNCTION CALL:
%MMN_child_FT(1.0,5) Hipass at 1 Hz and use the 5th standard
%
% OUTPUTS:

% THIS IS A WORK IN PROGRESS FUNCTION - any updates or suggestions would be
% much appreciated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MMN_child_FT(hipass) %hipass filter, which standard you want to use

downsample_factor = 5; %i.e. divide the orig sample rate by 5

global ft_default
ft_default.spmversion = 'spm12'; % Force SPM12, SPM8 doesn't go well with mac + 2017b
ft_defaults % This loads the rest of the defaults
ft_hastoolbox('fastica', 1);

files = dir('*.con');
hdr   = ft_read_header(files(1).name,'dataformat','yokogawa_con'); %hdr to get Fs etc.

if ~isempty(find(ismember(hdr.label,'AG160')))
    mkdir adultmeg
    fprintf('\n\nADULT MEG\n\n')
    error('This function is only for child MEG')
else
    fprintf('\n\nCHILD MEG\n\n')
end

% ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
cfg            = [];
cfg.trialfun   = 'ft_trialfun_general';
cfg.channel    = hdr.label(1:125);
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 40;
cfg.hpfilter   = 'yes';
cfg.hpfreq     = hipass; %FIXME should try something lower to appease reviewers
cfg.hpfiltord  = 5;
cfg.headerfile = files(1).name;
cfg.datafile   = files(1).name;
data           = ft_preprocessing(cfg);

cfg            = [];
cfg.resamplefs = data.fsample/downsample_factor; % Here we are downsampling from 1000Hz --> 200Hz
cfg.detrend    = 'yes'; % Helps with low-frequency drift
data           = ft_resampledata(cfg, data);

% deal with 50Hz line noise
cfg          = [];
cfg.bsfilter = 'yes';
cfg.bsfreq   = [49.5 50.5];
data         = ft_preprocessing(cfg, data);

data_clean = data;

%Create layout file for later + save
cfg      = [];
cfg.grad = data.grad; % struct containing gradiometer definition
lay      = ft_prepare_layout(cfg, data); % creates a 2-D layout of the channel locations

cfg                     = [];
cfg.dataset             = '2712_SH_ME125_2017_10_04_B1.con';
cfg.trialdef.prestim    = 1;         % pre-stimulus interval
cfg.trialdef.poststim   = 2;        % post-stimulus interval
cfg.trialdef.eventtype  = 'trigger';
cfg.trialdef.eventvalue = [1 2 3 4 5 6 7];%Trigger numbers
cfg.trialfun            = 'trig_fun_125_mmn';
cfg                     = ft_definetrial(cfg);

%Here we downsample the matrix representing the trial structure
cfg.trl(:,1) = round(cfg.trl(:,1)/downsample_factor); %FIXME this is possibly a bit dirty - is there a better way?
cfg.trl(:,2) = cfg.trl(:,1)+length(-(cfg.trialdef.prestim):1/data.fsample:cfg.trialdef.poststim)-1;
cfg.trl(:,3) = round(cfg.trl(:,3)/downsample_factor);
%*********************************************

%Epoch the data - first find trial indices corresponding to the events
%of interest and group
data          = ft_redefinetrial(cfg,data);
trigger_types = [cfg.event.value]';

cfg             = [];
cfg.method      = 'summary';
cfg.keepchannel = 'yes';
cfg.keeptrial   = 'nan';
data            = ft_rejectvisual(cfg, data);

deviant = find(ismember(trigger_types,1)); %trigger number 1 is first deviant

for i=1:5
    eval(['standard_',num2str(i),' = find(ismember(trigger_types,',num2str(i+1),'));']);
end

good_trials_idx = find(any(~isnan(cell2mat(cellfun(@(isgood)isgood(:),data.trial,'uni',0))),1)); % bad trial are NaN

for i=1:5
    cfg = [];
    eval(['cfg.trials   = intersect(standard_',num2str(i),',good_trials_idx);'])
    eval(['standard_',num2str(i),'_ave = ft_timelockanalysis(cfg,data);'])
end

cfg.trials  = intersect(deviant,good_trials_idx);
deviant_ave = ft_timelockanalysis(cfg,data);

c              = autumn;
cfg            = [];
cfg.layout     = lay;
cfg.xlim       = [-0.1 0.5];
cfg.graphcolor = ([c(1:10:50,:);[0 0 0]]);
figure;ft_multiplotER(cfg,standard_1_ave,standard_2_ave,standard_3_ave,standard_4_ave,standard_5_ave,deviant_ave);

cfg        = [];
cfg.method = 'power';
for i=1:5
    eval(['standard_',num2str(i),'_GFP = ft_globalmeanfield(cfg, standard_',num2str(i),'_ave);'])
end

deviant_GFP = ft_globalmeanfield(cfg, deviant_ave);

cfg       = [];
cfg.xlim  = [-0.1 0.5];
cfg.title = 'Global Field Power';
cfg.graphcolor = ([c(1:10:50,:);[0 0 0]]);
figure;subplot(4,2,1);ft_singleplotER(cfg,standard_1_GFP,standard_2_GFP,standard_3_GFP,standard_4_GFP,standard_5_GFP,deviant_GFP)
hold on
legend('standard_1','standard_2','standard_3','standard_4','standard_5','deviant')

cfg        = [];
cfg.method = 'triangulation';
neighbours = ft_prepare_neighbours(cfg, data.grad); %To get neigbours of the pak for plotting

%%find peak amplitude in M100 window and find the index - this will be used
%%for plotting
window_time = standard_5_GFP.time(standard_5_GFP.time>0.1&standard_5_GFP.time<0.2);
peak_time   = window_time(ismember(standard_5_GFP.avg(find(standard_5_GFP.time>0.1&standard_5_GFP.time<0.2)),...
    max(standard_5_GFP.avg(find(standard_5_GFP.time>0.1&standard_5_GFP.time<0.2)))));
peak_amp    = standard_5_GFP.avg(ismember(standard_5_GFP.avg(find(standard_5_GFP.time>0.1&standard_5_GFP.time<0.2)),...
    max(standard_5_GFP.avg(find(standard_5_GFP.time>0.1&standard_5_GFP.time<0.2)))));
peak_index  = find(ismember(standard_5_GFP.time,peak_time));
peak_amps   = standard_5_ave.avg(:,peak_index);
elec_idxs   = find(ismember(peak_amps,max(peak_amps)));


%%Plots below
cfg                  = [];
cfg.xlim             = [peak_time peak_time];
cfg.layout           = lay;
cfg.highlight        = 'on';
cfg.highlightchannel = vertcat(neighbours(elec_idxs).neighblabel,neighbours(elec_idxs).label);

subplot(4,2,[2 4]);ft_topoplotER(cfg,standard_5_ave)

cfg                = [];
cfg.channel        = vertcat(neighbours(elec_idxs).neighblabel,neighbours(elec_idxs).label);
cfg.avgoverchan    = 'yes';
cfg.nanmean        = 'yes';
standard_5_ave_ROI = ft_selectdata(cfg,standard_5_ave);
deviant_ave_ROI    = ft_selectdata(cfg,deviant_ave);

cfg       = [];
cfg.xlim  = [-0.1 0.5];
cfg.title = 'Peak Cluster';
subplot(4,2,3);ft_singleplotER(cfg,standard_5_ave_ROI,deviant_ave_ROI)
legend('standard','deviant')
end
