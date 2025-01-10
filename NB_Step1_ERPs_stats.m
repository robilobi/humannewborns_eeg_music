%% Analysis of erp to note onsets
% -- Roberta Bianco, Rome Nov 2022
% Load eeg, cut into mini epoch corresponding to the onset of each note
% Select events (25% highest and lowest IC).
% RUN IDyOM FIRST
% Plot the erp and run ttest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
U = 'C:\Users\robianco\OneDrive - Fondazione Istituto Italiano Tecnologia';
addpath([U '\MATools\fieldtrip-20220321\']); ft_defaults;
addpath('.\functions_RB');
addpath('.\templates');

%%% SET GRAPHICS / color palette
ft_hastoolbox('brewermap', 1);close all
cmap = colormap(flipud(brewermap(64,'RdBu')));
m = length(cmap);
cfg.ylim = [-0.55 0.55];
index1 = fix((0.4-cfg.ylim(1))/(cfg.ylim(2)-cfg.ylim(1))*m)+1; %A
index2 = fix((-0.4-cfg.ylim(1))/(cfg.ylim(2)-cfg.ylim(1))*m)+1; %A
% Convert indices to RGB values
RGB1 = cmap(index1, :);
RGB2 = cmap(index2, :);
RGB0 = [0.5, 0.5, 0.5];

RGB1 = ind2rgb(index1,cmap);
RGB2 = ind2rgb(index2,cmap);
RGB0 = RGB1; for i = 1:3; RGB0(:,:,i) = 0.5; end

%%% PATHS
stim_path = ([U '\BACHNB\STIMULI\']);
midi_path = ([U '\BACHNB\STIMULI\midiall\']);
midi_names = {'audio01.mid','audio02.mid','audio03.mid','audio04.mid','audio05.mid', ...
    'audio06.mid','audio07.mid','audio08.mid', 'audio09.mid','audio10.mid',...
    'shf01.mid','shf05.mid','shf08.mid','shf10.mid'};
code_corresponding_to_midi = [1:10 11 15 18 20];

%%% SETTINGS
dataINFolder  =  ([U '\BACHNB\RESULTS_NB\']);
features = {'Pitch' , 'Onset', 'ITI', 'IOI'};
for feat = 1:length(features)
    k = feat; 
    dataOUTFolder = ([dataINFolder 'ERP_1-30Hz_both+_' features{feat} '\']); mkdir(dataOUTFolder);

    %%% SUBJECT POOL
    eegsubjects = dir([dataINFolder,'*_PROC_1-30Hz*']);
    subjidx = 1:length(eegsubjects);
    computederp = 0;


    if computederp

        %% [all_id pitch onsets alliti allioi -So Eo - Sp Ep ];
        %      1     2     3      4       5     6  7     8  9 
        if      k == 1; predictor = 8; 
        elseif  k == 2; predictor = 6;
        elseif  k == 3; predictor = 4;
        elseif  k == 4; predictor = 5;
        end

        outname = 'FINAL_StimMatrix_both+';
        allinfo = readmatrix([stim_path 'allinfo_stimuli_' outname '.csv']);
        allinfo2 =[];
        midi_to_analyse = 1:14;
        for i = midi_to_analyse
            idx_midi_notes = ismember(allinfo(:,1), i);
            allinfo_thismel = allinfo(idx_midi_notes,:);

            svect = allinfo_thismel(:,predictor);   % S Sp and So values are the col 6:8
            Q    = quantile(svect,5);% IC 20% quantiles %use S overall 20 %for main anaysis, and Sp or So
            idx1 = (svect < Q(1));    % LOWEST IC
            allinfo_thismel(idx1,10) = 1; %10 is the column condition
            idx2 = (svect > Q(5));    % HIGHEST IC
            allinfo_thismel(idx2,10) = 2; %10 is the column condition
            allinfo_thismel(allinfo_thismel(:,4) == 99, 10) = 0;
            allinfo2 = [allinfo2; allinfo_thismel];
            quantile_20 = prctile(svect, 20);

        end

        %% EXTRACT MINIEPOCHS LOCKED TO ONSET PER CONDITION
        COND1 = cell(length(subjidx), 1);
        COND2 = cell(length(subjidx), 1);
        COND3 = cell(length(subjidx), 1);
        COND4 = cell(length(subjidx), 1);
        CONDR = cell(length(subjidx), 1);
        close all

        %% [all_id pitch onsets alliti allioi S0 E0 Sp Ep SpIn EpInt Spc Epc];
        %      1     2     3      4       5     6  7  8  9  10    11   12   13
        for sbj = 1:length(subjidx)
            file = [dataINFolder eegsubjects(sbj).name];
            [cond1, cond2, cond3, cond4, condrand, perc_trremoved]  = fun_load_epoch_eeg_erp_new(file, allinfo2, midi_names,code_corresponding_to_midi, sbj);
            COND1{sbj}=cond1; %L or
            COND2{sbj}=cond2; %H or
            COND3{sbj}=cond3; %L sh
            COND4{sbj}=cond4; %H sh
            CONDR{sbj}=condrand; % random trigger
            clear cond1 cond2 cond3 cond4 condrand
            trialremoved(sbj,:) = perc_trremoved;
        end

        save([dataOUTFolder 'COND1'], 'COND1'); %Real Low
        save([dataOUTFolder 'COND2'], 'COND2'); %Real High
        save([dataOUTFolder 'COND3'], 'COND3'); %Shuf Low
        save([dataOUTFolder 'COND4'], 'COND4'); %Shuf High
        save([dataOUTFolder 'CONDR'], 'CONDR'); %
        save([dataOUTFolder 'trialremoved_persubj'], 'trialremoved'); %

    else

        load([dataOUTFolder 'COND1']); %Real Low
        load([dataOUTFolder 'COND2']); %Real High
        load([dataOUTFolder 'COND3']); %Shuf Low
        load([dataOUTFolder 'COND4']); %Shuf High
        load([dataOUTFolder 'CONDR']); %
    end


    subjidx = [1:10 12:33 35:44 46:47 54:58];

    %% TRIM DATA
    count = 1;
    for sbj = subjidx
        cfg = [];
        cfg.latency = [-0.05 0.400];
        trim_COND1{count}=ft_selectdata(cfg,COND1{sbj});
        trim_COND2{count}=ft_selectdata(cfg,COND2{sbj});
        trim_COND3{count}=ft_selectdata(cfg,COND3{sbj});
        trim_COND4{count}=ft_selectdata(cfg,COND4{sbj});
        trim_CONDR{count}=ft_selectdata(cfg,CONDR{sbj});
        count = count+1;
    end


    %% GRANDAVERAGE

    cfg = [];
    cfg.channel = 'all';
    cfg.latency = 'all';
    cfg.parameter = 'avg';
    grandavg_cond1 = ft_timelockgrandaverage(cfg, trim_COND1{:});
    grandavg_cond2 = ft_timelockgrandaverage(cfg, trim_COND2{:});
    grandavg_cond3 = ft_timelockgrandaverage(cfg, trim_COND3{:});
    grandavg_cond4 = ft_timelockgrandaverage(cfg, trim_COND4{:});
    grandavg_condr = ft_timelockgrandaverage(cfg, trim_CONDR{:});


    %% CLUSTER BASED STATISTICS
    load('Layout_Newborn_EEG.mat');
    cfg = [];
    cfg.layout          = lay;
    cfg.skipscale       = 'yes'; cfg.skipcomnt   = 'yes' ;
    layout              = ft_prepare_layout(cfg,grandavg_cond1);
    cfg                 = [];
    cfg.layout          = layout;
    cfg.method          = 'distance'; % for prepare_neigh
    cfg.neighbourdist   = 0.15;         % results in avg 5 channels
    neighbours          = ft_prepare_neighbours(cfg);
    cfg = [];
    cfg.channel          = grandavg_cond1.label;
    cfg.method           = 'ft_statistics_montecarlo';  % use the Monte Carlo method to calculate probabilities
    cfg.statistic        = 'ft_statfun_depsamplesT';    % use the dependent samples T-statistic as a measure to evaluate the effect at each sample
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.1;                        % threshold for the sample-specific test, is used for thresholding
    cfg.clusterstatistic = 'maxsum';
    cfg.clusterthreshold = 'nonparametric_common';
    cfg.minnbchan        = 3;                           % minimum number of neighbouring channels that is required
    cfg.tail             = 0;                           % test the left, right or both tails of the distribution
    cfg.clustertail      = cfg.tail;
    cfg.alpha            = 0.05;                        % alpha level of the permutation test
    cfg.correcttail      = 'alpha';                     % see https://www.fieldtriptoolbox.org/faq/why_should_i_use_the_cfg.correcttail_option_when_using_statistics_montecarlo/
    cfg.computeprob      = 'yes';
    cfg.numrandomization = 1000;                        % number of random permutations
    cfg.neighbours       = neighbours;                  % the neighbours for each sensor to form clusters
    Nsub = length(subjidx);
    cfg.design(1, 1:2*Nsub) = [ones(1, Nsub) 2*ones(1, Nsub)];
    cfg.design(2, 1:2*Nsub) = [1:Nsub 1:Nsub];
    cfg.ivar = 1; % the 1st row in cfg.design contains the independent variable
    cfg.uvar = 2; % the 2nd row in cfg.design contains the subject number
    stat_cluster_re = ft_timelockstatistics(cfg, trim_COND2{:}, trim_COND1{:});
    filename=([dataOUTFolder 'ClstPerm_image_RE']);
    fun_clusterBasedPermutation_Plot(stat_cluster_re,filename, 'C2vsC1', 'images','F_value', layout);
    stat_cluster_sh = ft_timelockstatistics(cfg, trim_COND4{:}, trim_COND3{:});
    filename=([dataOUTFolder 'ClstPerm_image_SH']);
    fun_clusterBasedPermutation_Plot(stat_cluster_sh,filename, 'C2vsC1', 'images','F_value', layout);

    %% plots clusters

    lay_mk= load('Layout_Monkey_EEG.mat');
    order_chan = lay_mk.lay.label;

    h=figure;clf
    h.Position = [100 100 800 700];
    condstat(1) = stat_cluster_re;
    condstat(2) = stat_cluster_sh;
    titles = {'Real', 'Shuffled'};
    for i = 1:2
        stat = condstat(i);
        time = stat.time;
        m1_label = stat.label;
        [a, ord] = ismember(order_chan, m1_label);
        ord(ord==0) =[];
        m1 = stat.mask(ord,:);

        limits = [-6 6];
        if ~isfield(stat, 'posclusterslabelmat')
            stat.posclusterslabelmat = zeros(size(ord,1),length(time));end
        if ~isfield(stat, 'negclusterslabelmat')
            stat.negclusterslabelmat = zeros(size(order_chan,1),length(time));end
        M = stat.mask(ord,:).*(stat.posclusterslabelmat(ord,:)-stat.negclusterslabelmat(ord,:));
        subplot(1,2,i);
        imagesc(time,1:length(order_chan),M, limits);hold on;
        yticks(1:size(order_chan,1));
        yticklabels(order_chan); title(titles{i})
        colormap(flipud(brewermap(64,'RdBu')));

    end
    filename=([dataOUTFolder 'ClstPerm_2D_realvsshf']);
    saveas(gcf, [filename '.png']);
    saveas(gcf, [filename '.fig']);
    print([filename '.pdf'],'-dpdf','-bestfit')

    %% Plot ERP and Topo
    climd = min(min([mean(grandavg_cond1.avg)]))-0.08;
    climu = max(max([mean(grandavg_cond3.avg)]))+0.08;

    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    Oc1 = ft_math(cfg,grandavg_cond1,grandavg_condr);
    Oc2 = ft_math(cfg,grandavg_cond2,grandavg_condr);
    Sc1 = ft_math(cfg,grandavg_cond3,grandavg_condr);
    Sc2 = ft_math(cfg,grandavg_cond4,grandavg_condr);

    REdif = ft_math(cfg,grandavg_cond2,grandavg_cond1);
    SHdif = ft_math(cfg,grandavg_cond4,grandavg_cond3);


    %% Non fieldtrip plotting
    %%% pull data in non fieldtrip structure for better plotting
    allData = [];
    if feat == 1 || feat == 2;  chantoplot = {'AF3', 'AF4', 'AFz', ...
            'F1', 'Fz', 'F2', ...
            'FC1', 'FC2', 'FC3', 'FC4', ...
            'Fp1', 'Fp2', 'AF7', 'AF8'};
    elseif feat == 4 || feat == 3; chantoplot = { 'F1', 'Fz', 'F2', ...
            'FC1', 'FC2', 'FC3', 'FC4', ...
            'C2', 'C1', 'C3', 'C4',...
            'CP1', 'CPz', 'CP2'};
    end
    myString = strjoin(cellfun(@string, chantoplot), ',');
    conditions = {'Low (real)', 'High (real)', 'Low (sh)', 'High (sh)'};
    for si = 1:length(subjidx)
        [a,ch] = ismember( chantoplot, trim_COND1{si}.label);
        allData(si,1,:,:) = trim_COND1{si}.avg(ch, :); % sub x cond x chn x time
        allData(si,2,:,:) = trim_COND2{si}.avg(ch, :);
        allData(si,3,:,:) = trim_COND3{si}.avg(ch, :);
        allData(si,4,:,:) = trim_COND4{si}.avg(ch, :);
    end
    B = 500;
    bsmean =[]; bsstd=[];stderror=[];
    for c = 1:4
        x = (allData(:,c,:,:));
        x = squeeze(mean(x,3));
        %  Force x to be time*repetitions
        x=x';
        [mn,sd,bsall] = fBootstrapMean(x,B);
        bsmean(c,:)=mn';
        bsstd(c,:)=sd'; % Bootstrap standard deviation IS the estimate of the sample standard ERROR!
        stderror(c,:) = std(x,0,2)/sqrt(numel(subjidx));
    end

    databoot = squeeze(mean(allData, 3));
    diffRE=squeeze(databoot(:,2,:)-databoot(:,1,:)); % sub x time
    diffSH=squeeze(databoot(:,4,:)-databoot(:,3,:));
    dataB=bootstrap(diffRE);
    sRE=findSigDiff(dataB, 0.05);
    dataB=bootstrap(diffSH);
    sSH=findSigDiff(dataB, 0.05);

    RGB1 = cmap(index1, :); RGB2 = cmap(index2, :);
    C=[RGB2;
        RGB1;
        RGB2;
        RGB1];
    time_s = grandavg_cond1.time;

    %%%% Plot boostrap SD and mean across subj for all conditions
    fh=figure(1);clf
    fh.Position = [0, 0, 600, 600];
    subplot(2,4, 1:2);
    for c = 1:2; plot(time_s,bsmean(c,:),'color',C(c,:), 'LineWidth', 1.5); hold on; end
    xlim([-0.05 0.4]);
    % Plot SE of the bootstrap resampling
    for c = 1:2 % loop over conditions
        b = bsstd(c,:)'; % STDEV
        a = bsmean(c,:)';
        Y = [b+a;  flipud(-b+a)]';
        abscissa = time_s';
        X = [abscissa; flipud(abscissa)];
        h = fill(X,Y,C(c,:),'edgecolor','none','facealpha',0.2); hold on;
        plot(abscissa,a*0,'-k'); % plot zero line
    end
    ylim([-0.2 0.2]);
    plot(time_s, -0.08*abs(sRE),'k', LineWidth=3);
    lh= legend(conditions(1:2), 'Location','southeast');
    set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
    title(['\bf p=.01 uncorr, REAL, chn:' myString{1}]); % \bf = bold font
    % title(['\bf p=.05, chn: ' myString{1}]); % \bf = bold font
    set(gca,'Box','off'); ylabel('\bf Voltage (uV)');

    subplot(2,4,3:4);
    for c = 3:4; plot(time_s,bsmean(c,:),'color',C(c,:), 'LineWidth', 1.5); hold on; end
    xlim([-0.05 0.4]);
    % Plot SE of the bootstrap resampling
    for c = 3:4 % loop over conditions
        b = bsstd(c,:)'; % STDEV
        a = bsmean(c,:)';
        Y = [b+a;  flipud(-b+a)]';
        abscissa = time_s';
        X = [abscissa; flipud(abscissa)];
        h = fill(X,Y,C(c,:),'edgecolor','none','facealpha',0.2); hold on;
        plot(abscissa,a*0,'-k'); % plot zero line
    end
    ylim([-0.2 0.2]);
    plot(time_s, -0.08*abs(sSH),'k', LineWidth=3);
    lh= legend(conditions(3:4), 'Location','southeast');
    set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
    %title(['\bf p=.01 uncorr, SHUF, chn:' myString{1}]); % \bf = bold font
    set(gca,'Box','off'); ylabel('\bf Voltage (uV)');


    cfg = [];
    cfg.layout = layout;
    cfg.figure = 'gca';
    cfg.colorbar = 'no';
    cfg.commentpos = 'middletop';
    cfg.marker = 'on';
    cfg.style = 'both'; % 'straight';
    cfg.ylim = [climd climu];
    cfg.zlim = [-0.1 0.1];
    cfg.xlim = [0.16 0.18];
    cfg.comment = num2str(cfg.xlim);
    subplot(2,4,5); ft_topoplotER(cfg,REdif);title([' H vs L (real)']);
    subplot(2,4,7); ft_topoplotER(cfg,SHdif);title([' H vs L (sh)']);
    cfg.xlim = [0.31 0.36];
    cfg.comment = num2str(cfg.xlim);cfg.colorbar = 'no';
    subplot(2,4,6); ft_topoplotER(cfg,REdif);title([' H vs L (real)']);
    subplot(2,4,8); ft_topoplotER(cfg,SHdif);title([' H vs L (sh)']);
    c = colorbar;c.LineWidth = 0.1;
    c.Position = [0.92, 0.1, 0.01, 0.08]; % [x, y, width, height]
    c.FontSize = 8;
    colormap(flipud(brewermap(64,'RdBu')));

    filename=([dataOUTFolder 'FIG_ERP_TOPO']);
    saveas(gcf, [filename '.png']);
    print([filename '.pdf'],'-dpdf','-bestfit')
    print([filename '.eps'],'-depsc')

end