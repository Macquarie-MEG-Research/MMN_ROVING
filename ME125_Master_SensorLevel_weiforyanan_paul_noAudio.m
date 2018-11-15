%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This script is created for analyzing %%%%%%
%%%%%% Roving MMN data acquaried at KIT-MQ  %%%%%%
%%%%%%    funded by ARC-DP [DP170103148]    %%%%%%
%%%%%%  prepared by Wei He using fieldtrip  %%%%%%
%%%%%% and functions written by Paul Sowman %%%%%%
%%%%%% and Robert Seymour , August 2018     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%!!! eash subject folder should have 5 files: .con .hsp .elp .mrk .txt(matlab trigger output) !!!
%%!!! manually reformat the .mrk name to '*_ini.mrk' !!!
%%!!! make sure Fieldtrip toolbox is in the path !!!
%%!!! make sure "scripts" folder is in the path !!!
%%!!! make sure "database_for_MEMES" folder is in the path !!!


%%
clc
close all

% % add fieldtrip in path
% restoredefaultpath
% addpath E:\fieldtrip-20180825; % change path if necessary
% ft_defaults
% fprintf('\nAdd Fieldtrip into path.\n')
%
% % add scripts full path
% addpath(genpath('E:\ARC-DP\RawData\MMN\scripts'));% change path if necessary
% fprintf('\nAdd scripts into path.\n')
%
% %% BATCH STARTS
%
% subjects = input('please type the number of subjects\n','s');
% ee = str2num(subjects);%offset number of directories
% y=subdir(cd);

y = dir('*rethm.con');

for i = 1:length(y)
    x = [];
    close all;
    if strcmp(y.name(end-8:end-4), 'rethm')
        %cd(y{i});
        %fprintf('\n%s is a ReTHM directory! Perform analysis\n', y{i})
        
        %% PART 1: MEMES:  Co-registration between MEG, Head Surface, and Brain Sruface in the MRI-library
        fname    = dir('*.elp');
        fname    = fname.name(1:strfind(fname.name,'.')-1);
        pathname = [pwd,'/'];
        elpfile  = [fname,'.elp'];
        hspfile  = [fname,'.hsp'];
        confile  = [fname,'_B1_denoise_rethm.con']; % !!! if NOT ReTHM data change the con file name here !!!
        mrkfile  = [fname,'_INI.mrk'];
        bad_coil = [];
        %path_to_MRI_library         = 'C:\Users\mq42604613\Desktop\ARC_DP_pipeline_test\database_for_MEMES\'; % !!! manually check !!!
        
        %         child_MEMES(pathname,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,bad_coil,'yes');  % found here that
        
        %% PART 2: Functional processing - Sensor Level ERFs - pre-&post-Deviants
        hdr      = ft_read_header(confile,'dataformat','yokogawa_con'); %hdr to get Fs etc.
        filename = confile;
        
        %%2.1 Raw data & Bandpass Filter
        cfg            = [];
        cfg.trialfun   = 'ft_trialfun_general';  %ft_definetrial: defines the segments of data that will be read in by FT_PREPROCESSING
        cfg.channel    = hdr.label(1:125);
        cfg.lpfilter   = 'yes';
        cfg.lpfreq     = 40;
        cfg.hpfilter   = 'yes';
        cfg.hpfreq     = 1; %FIXME should try something lower to appease reviewers
        cfg.hpfiltord  = 5;
        cfg.dftfreq    = 50; % removal line noise
        cfg.headerfile = filename;
        cfg.datafile   = filename;
        data           = ft_preprocessing(cfg);
        
        %%2.2 Create layout file for later + save
        cfg      = [];
        cfg.grad = data.grad; % struct containing gradiometer definition
        lay      = ft_prepare_layout(cfg, data); % creates a 2-D layout of the channel locations
        %ft_layoutplot(cfg);
        save (['lay_', num2str(fname(1:4)),'.mat'],'lay')
        
        %%2.3 Trigger-based trial selection
        cfg                     = [];
        cfg.dataset             = filename;
        cfg.path                = pathname;
        cfg.trialdef.prestim    = 0.1;         % pre-stimulus interval
        cfg.trialdef.poststim   = 0.5;        % post-stimulus interval
        cfg.trialdef.eventtype  = 'trigger';
        cfg.trialdef.eventvalue = 1:7;%Trigger numbers
        cfg.Fs                  = hdr.Fs;
        cfg.First_Channel       = 146; % First trigger channel
        cfg.Last_Channel        = 158; % Last trigger channel
        cfg.Audio_Channel       = 0; % 0 for none
        cfg.fixed_offset        = 42;
        cfg.trialfun            = 'FindTriggers_AudioCh';
        cfg                     = ft_definetrial(cfg);
        
        alldata = ft_redefinetrial(cfg,data);
        event   = cfg.event;
        trl     = cfg.trl;
        
        save (['event_', num2str(fname(1:4)),'.mat'], 'event')
        save (['trl_', num2str(fname(1:4)),'.mat'], 'trl')
        
        % Detrend and demean each trial
        cfg         = [];
        cfg.demean  = 'yes';
        cfg.detrend = 'yes';
        % cfg.baselinewindow          = [-0.5 0];
        alldata     = ft_preprocessing(cfg,alldata);
        
        %%2.4 Visual artefact rejection
        [z,bad_trials,data_clean] = artifacts_max_z(alldata,10);
        
        good_trials_idx = setdiff(1:length(event),bad_trials);
        
        cfg        = [];
        cfg.trials = good_trials_idx;
        data_clean = ft_selectdata(cfg,data_clean); %remove bad trials
        
        event_clean=event(good_trials_idx); %remove bad trials events
        trl_clean=trl(good_trials_idx,:);
        
        %good_trials_idx             = find(~isnan(cell2mat(cellfun(@(isgood)isgood(1),data_clean.trial,'uni',0)))); %just need to evaluate the first element as all samples in bad trial are NaN
        %         for j = 1:length(bad_trials)
        %             event(bad_trials(j)).bad = 1;
        %         end
        event=event_clean;
        
        for j=1:length(event)
            if event(j).value == 1
                event(j).type = 'deviant';
            else
                event(j).type = 'standard';
            end
        end
        
        for j=1:length(event)-1
            if strcmp(event(j).type,'standard') & strcmp(event(j+1).type,'deviant')
                event(j).type   = 'predeviant';
                event(j+1).type = event(j+1).type;
            else
            end
        end
        
        %%2.5 Downsample
        downsample_factor = 5;
        
        cfg            = [];
        cfg.resamplefs = data_clean.fsample/downsample_factor; % Here we are downsampling from 1000Hz --> 200Hz: Garrido et al., 2008, Neuroimage
        cfg.detrend    = 'yes'; % Helps with low-frequency drift
        data_clean     = ft_resampledata(cfg, data_clean);
        
        save (['data_clean_', num2str(fname(1:4)),'.mat'], 'data_clean')
        
        %%2.6 Epoching
        trigger_types  = {event.type};
        trigger_values = cell2mat({event.value});
        
        % pre-deviants/deviants
        predeviant_trials = find(ismember(trigger_types,'predeviant')); %trigger number 1 is first deviant
        deviant_trials    = find(ismember(trigger_types,'deviant')); %trigger number 1 is first deviant
        
        %         repeat_1 = setdiff(find(ismember(trigger_values,2)),bad_trials); %trigger number
        %         repeat_2 = setdiff(find(ismember(trigger_values,3)),bad_trials); %trigger number
        %         repeat_3 = setdiff(find(ismember(trigger_values,4)),bad_trials); %trigger number
        %         repeat_4 = setdiff(find(ismember(trigger_values,5)),bad_trials); %trigger number
        %         repeat_5 = setdiff(find(ismember(trigger_values,6)),bad_trials); %trigger number
        %         repeat_6 = setdiff(find(ismember(trigger_values,7)),bad_trials); %trigger number
        %
        %         %equalise repeat numbers
        %         repeat_1 = repeat_1(randperm(length(repeat_1)));
        %         repeat_1 = repeat_1(1:length(repeat_6));
        %
        %         repeat_2 = repeat_2(randperm(length(repeat_2)));
        %         repeat_2 = repeat_2(1:length(repeat_6));
        %
        %         repeat_3 = repeat_3(randperm(length(repeat_3)));
        %         repeat_3 = repeat_3(1:length(repeat_6));
        %
        %         repeat_4 = repeat_4(randperm(length(repeat_4)));
        %         repeat_4 = repeat_4(1:length(repeat_6));
        %
        %         repeat_5 = repeat_5(randperm(length(repeat_5)));
        %         repeat_5 = repeat_5(1:length(repeat_6));
        
        cfg        = [];
        cfg.trials = setdiff(deviant_trials,bad_trials);
        deviant    = ft_redefinetrial(cfg,data_clean);
        save (['deviant_', num2str(fname(1:4)),'.mat'], 'deviant')
        
        cfg        = [];
        cfg.trials = setdiff(predeviant_trials,bad_trials);
        predeviant = ft_redefinetrial(cfg,data_clean);
        save (['predeviant_', num2str(fname(1:4)),'.mat'], 'predeviant')
        
        
        %%2.7 Averaging
        cfg            = [];
        deviant_ave    = ft_timelockanalysis(cfg,deviant);
        predeviant_ave = ft_timelockanalysis(cfg,predeviant);
        save (['deviant_ave_', num2str(fname(1:4)),'.mat'], 'deviant_ave')
        save (['predeviant_ave_', num2str(fname(1:4)),'.mat'], 'predeviant_ave')
        
        
        %%2.8 Planar gradient transform
        % Calculate the planar gradient of the averaged data:
        cfg            = [];
        cfg.method     = 'triangulation';
        cfg.neighbours = ft_prepare_neighbours(cfg, data_clean.grad); %To get neigbours of the pak for plotting
        
        cfg.planarmethod      = 'sincos';
        deviant_ave_planar    = ft_megplanar(cfg, deviant_ave);
        predeviant_ave_planar = ft_megplanar(cfg, predeviant_ave);
        
        
        % Combine the horizontal and vertical components of the planar gradient
        cfg                        = [];
        deviant_ave_planar_comb    = ft_combineplanar(cfg,deviant_ave_planar);
        cfg                        = [];
        predeviant_ave_planar_comb = ft_combineplanar(cfg,predeviant_ave_planar);
        save (['deviant_ave_planar_comb_', num2str(fname(1:4)),'.mat'], 'deviant_ave_planar_comb')
        save (['predeviant_ave_planar_comb_', num2str(fname(1:4)),'.mat'], 'predeviant_ave_planar_comb')
        
        
        % GFP for planar gradient
        cfg                   = [];
        cfg.method            = 'power';
        deviant_planar_GFP    = ft_globalmeanfield(cfg, deviant_ave_planar_comb);
        predeviant_planar_GFP = ft_globalmeanfield(cfg, predeviant_ave_planar_comb);
        
        cfg                       = [];
        cfg.operation             = 'subtract';
        cfg.parameter             = 'avg';
        difference_planar_GFP     = ft_math(cfg, deviant_planar_GFP , predeviant_planar_GFP );
        difference_planar_GFP.avg = abs(difference_planar_GFP.avg);
        %
        %         cfg            = [];
        %         cfg.xlim       = [-0.1 0.4];
        %         cfg.title      = 'Global Field Power';
        %         cfg.graphcolor = 'brk';
        %         figure;
        %         % ft_singleplotER(cfg,standard_1_GFP,standard_2_GFP,standard_3_GFP,standard_4_GFP,standard_5_GFP,deviant_GFP)
        %         ft_singleplotER(cfg,predeviant_planar_GFP ,deviant_planar_GFP , difference_planar_GFP)
        %         legend('pre-deviant','deviant', 'difference')
        %         title(['GPF of the difference waveform in subject ', num2str(fname(1:4)),' (planar gradients)']); drawnow;
        %         print(['GFP_StdvsDev_', num2str(fname(1:4))],'-dpng');
        %
        %         % Toppo plots on the deviant waveform (planar gradients)
        %         % To plot a squence of topo plots equally spaced between 0.1 and 0.4 second
        %         timestep      = 0.05; % in 50 ms steps
        %         sampling_rate = data_clean.fsample;
        %         sample_count  = length(deviant_ave_planar_comb.time);
        %
        %         j = [0:timestep:0.4];
        %         % m                           = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
        %         ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
        %
        %         figure;
        %         for k = 1: length(j)-1
        %             cfg         = [];
        %             cfg.comment = 'no';
        %             %     cfg.marker          = 'off';
        %             cfg.layout  = lay;
        %             %     cfg.colorbar        = 'southoutside';
        %             cfg.style   = 'straight';
        %
        %             subplot(2,round(length(j)-1)/2,k)
        %             cfg.xlim = [j(k) j(k+1)];
        %             cfg.zlim = 'maxabs';
        %             ft_topoplotER(cfg, deviant_ave_planar_comb)
        %             colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
        %             %     colorbar
        %             title(['Time [', num2str(j(k)),' ', num2str(j(k+1)),']']); drawnow;
        %             hold on
        %         end
        %         h = suptitle (['Topographic plots of the deviant response in subject ', num2str(fname(1:4)),' (planar gradients)']);
        %         set (h,'FontSize',12,'FontWeight','bold')
        %         set(gcf, 'Position', [200, 200, 1000, 1000])
        %         print(['Topo_deviantavg_planar_', num2str(fname(1:4))],'-dpng');
        %
        %         figure;
        %         for k = 1: length(j)-1
        %             cfg         = [];
        %             cfg.comment = 'no';
        %             %     cfg.marker          = 'off';
        %             cfg.layout  = lay;
        %             %     cfg.colorbar        = 'southoutside';
        %             cfg.style   = 'straight';
        %
        %             subplot(2,round(length(j)-1)/2,k)
        %             cfg.xlim = [j(k) j(k+1)];
        %             cfg.zlim = 'maxabs';
        %             ft_topoplotER(cfg, predeviant_ave_planar_comb)
        %             colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
        %             %     colorbar
        %             title(['Time [', num2str(j(k)),' ', num2str(j(k+1)),']']); drawnow;
        %             hold on
        %         end
        %         h = suptitle (['Topographic plots of the standard response in subject ', num2str(fname(1:4)),' (planar gradients)']);
        %         set (h,'FontSize',12,'FontWeight','bold')
        %         set(gcf, 'Position', [200, 200, 1000, 1000])
        %         print(['Topo_standardavg_planar_', num2str(fname(1:4))],'-dpng');
        
    else
        %        fprintf('\n%s is NOT a ReTHM directory!\n', y{i})
    end
end















