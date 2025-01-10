%% PLot TRF output%
clear all;
close all;clc;
clear all;
U = 'C:\Users\robianco\OneDrive - Fondazione Istituto Italiano Tecnologia';
addpath([U '\MATools\fieldtrip-20220321\']); ft_defaults;
addpath('.\functions_RB')
addpath('.\templates')
ft_defaults;
ft_colormap('RdBu');


dataOUTfolder = [U '\BACHNB\RESULTS_NB\TRF_GT_both+_rnd_sosp\'];
models={'full','-Mus', '-So', '-Sp', '-IOI', 'ITI'};

%%% Adjust layout (remove electrodes)
cfg = [];
cfg.layout = 'Layout_Newborn_EEG.mat';
cfg.baselinetype = 'absolute';
cfg.skipscale   = 'yes'; cfg.skipcomnt   = 'yes' ;
layout = ft_prepare_layout(cfg);
norm_chanlbls =layout.label([2:26 28:65]);
sel = ismember(layout.label, norm_chanlbls);
layout.pos = layout.pos(sel,:);
layout.height = layout.height(sel);
layout.label = layout.label(sel);
layout.width = layout.width(sel);

%%%% STIMULUS names
stimnames = {'audio01','audio02','audio03','audio04','audio05', ...
    'audio06','audio07','audio08', 'audio09','audio10',...
    'shf01','shf05','shf08','shf10'};

%%%% LOAD DATA
fn= 'TRFout_newborns_AM.mat';
load([dataOUTfolder fn]);
nmodel = length(model);
subjidx = [1:10 12:33 35:44 46:47 54:58];

ROI(1).chn= {}; %
% ROI(1).chn= {'Fp1', 'AF7','AF3', 'F9', 'F7', 'F5', 'F3' , 'FT7', 'FC5','FC3', 'T7', 'C5' , 'C3'}; %Anterior left
% ROI(2).chn= {'Fp2', 'AF8','AF4', 'F10','F8', 'F6', 'F4' , 'FT8', 'FC6','FC4', 'T8', 'C6' , 'C4'}; %Anterior right
% ROI(3).chn= {'AFz', 'Fz','F1', 'F2', 'FC1', 'FC2' , 'C1', 'C2','CPz', 'CP1', 'CP2'};      %Anterior central
% ROI(4).chn= {'TP7', 'CP5','CP3', 'CP1', 'CP4', 'CP6', 'TP8' , 'P9', 'P7','P5', 'P3', 'P1', 'Pz',...
%     'P2', 'P4', 'P6', 'P8', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'O1' , 'Oz', 'O2'};      %Posterior
% for k = 1:length(ROI)
%     [a,ch] = ismember( ROI(k).chn, norm_chanlbls);
%     idxroi(k).idxch = ch;
% end

Nsbj = length(subjidx);
Matrixfinal = [];


%% Select the channels (best 25% in the full model) and subjects
exclude = [];
for sbj   = subjidx
    mm1 = mean(model(1).statsall{sbj}.r_test);
    mm = mm1;
    [sorted, idx]= sort(mm, 'descend');    
    idxch = idx(sorted>0);
    perchntokeep = 25;
    chntokeep = ceil((63*perchntokeep)/100);
    if length(idxch) < chntokeep
        badsubj = sbj;
        exclude = [exclude badsubj];
    end    
    if length(idxch) > chntokeep; idxch = idxch(1:chntokeep);end
%     if length(idxch) == 0 ; idxch = idx(1);end
    idxroi(1).idxch = idxch;
    idxselectedchannles(sbj).idxch = idxch;
    idxselectedchannles(sbj).count = length(idxch);

end
save([dataOUTfolder 'idxselectedchannles.mat'], 'idxselectedchannles', '-mat');
%subjidx = setdiff(subjidx, exclude);

% figure(10); set(gcf,'Position',[100 100 1600 800]);
% count =1;
% for sbj   = subjidx
%     mm1 = mean(model(1).statsall{sbj}.r_test);
%     mm2 = mean(model(2).statsall{sbj}.r_test);
%     d = mm1-mm2;
%     d = mm1;
%     climd = -abs(max([nanmedian(d)]));
%     climu = abs(max([nanmedian(d)]));
%     subplot(7,7,count); ft_plot_topo(layout.pos(:,1),layout.pos(:,2),d ,...
%         'mask',layout.mask,'outline',layout.outline, 'interplim','mask','clim',[climd climu]);
%     set(gca,'FontSize',9);axis off; axis('square'); colorbar;
%     colormap(flipud(brewermap(64,'RdBu')));title(num2str(sbj));
%     count = count +1;
% end
% sgtitle('AM model');
% fn = sprintf(['/FIG_GAVG_Rtopo_ALL_AMbysbj_' ]);
% saveas(gcf,[dataOUTfolder fn '.png']);
% print([dataOUTfolder fn '.pdf'],'-dpdf','-bestfit')

%% Start
count = 1;
for sbj   = subjidx
    % some sbj don't have all
    songs_played = find(ismember(stimnames, model(1).statsall{1,sbj}.stim_names ));
    idxch = idxselectedchannles(sbj).idxch;

    % select channels: sbj x melody x roi (mean over channels)
    for i = 1:nmodel
        for j = 1:length(ROI)
            %             idxch = idxroi(j).idxch;
            r(i).model(count, songs_played,j) = mean(model(i).statsall{sbj}.r_test(:,idxch), 2) ;
        end
    end

    % grandaverage of r values subj x melody x channels
    for i = 1:nmodel
        gavgr(i).model(count, songs_played, :)   = model(i).statsall{sbj}.r_test;
    end

    % Store for grandaverage of weights subj x pred x time x channels
    for i=1:nmodel
        avgmodel(i).avg = mTRFmodelAvgRB(model(i).mdl{sbj},1);
        gavgmodel(i).avg(count, :, :, :)   = avgmodel(i).avg.w;
    end

    % grandaverage of weights subj x melody x pred x time x channels
    for i=1:nmodel
        gavgw(i).model(count, songs_played, :, :,:) = cellfun(@(x) x.w, model(i).mdl{sbj}, 'UniformOutput', false);
    end


    % EXPORT TABLE FOR STATS IN r
    for j = 1:length(ROI)
        idxch = idxroi(j).idxch;
        nsong= length(songs_played);
        AO = r(1).model(count,songs_played, j);
        A = [AO, AO,AO,AO, AO];
        iti = AO - r(6).model(count,songs_played, j);
        ioi = AO - r(5).model(count,songs_played, j);
        sp = AO - r(4).model(count,songs_played, j);
        so = AO - r(3).model(count,songs_played, j);
        mus = AO - r(2).model(count,songs_played, j);

        n = nmodel-1;
        names = stimnames(1:nsong);
        rmat = [so, sp, ioi,iti, mus];

        cond = repmat([1 1 1 1 1 1 1 1 1 1 2 2 2 2], 1, n);
        modelname = [];roinumber=[];
        for k = 1:n
            modelname = [modelname, repmat(k,1,nsong) ];
        end
        [g, h]=ismember(names, stimnames);
        stimid =repmat(h, 1,n) ;
        melbymodel = nsong*n;
        roinumber = repmat(j, 1,melbymodel);
        Matrix(1:melbymodel,1) = sbj;                   %subj
        Matrix(1:melbymodel,2) = rmat;                  %Pred acc difference Am / Amc - A
        Matrix(1:melbymodel,3) = cond;                  %Or or Sh Mel
        Matrix(1:melbymodel,4) = modelname;             %model AM, AMp, AMo
        Matrix(1:melbymodel,5) = stimid;                %stim ID
        Matrix(1:melbymodel,6) = A;                     %Pred Acc Full model
        Matrix(1:melbymodel,7) = roinumber;             %N ROI
        Matrixfinal = [Matrixfinal; Matrix];
        Matrix = [];
    end

    count = count +1;
    clear A
end
csvwrite([dataOUTfolder 'PredAccTRFforstat.csv'], Matrixfinal)


%% PLOT BARS Individual subjects
all = [];
k = 1; %ROI
goodsbj = ismember(subjidx, subjidx);
figure(97); set(gcf,'Position',[100 100 1500 800]);
all = mean(r(1).model(:,:, k) ,2);
[a,b]=sort(all, 'descend');
subplot(3,nmodel-1, 1:nmodel-1);
bar(all(b)); title('AM');xlabel('subjects');
xticks(1:length(b)); xticklabels(num2str(subjidx(b)'));
p=nmodel;
for i=2:nmodel
    subplot(3,nmodel-1,p);
    real = mean(r(1).model(goodsbj,1:10, k) ,2) - mean(r(i).model(goodsbj,1:10, k) ,2);
    bar(real(b));title(sprintf(models{i})); ylim([-0.009 0.01]);
    subplot(3,nmodel-1,p+nmodel-1);
    real = mean(r(1).model(goodsbj,11:14, k) ,2) - mean(r(i).model(goodsbj,11:14, k) ,2);
    bar(real(b));title(sprintf(models{i})); ylim([-0.009 0.01]);
    p = p+1;
end
fn = sprintf(['/FIG_BAR_fullmodel_' ]);
saveas(gcf,[dataOUTfolder fn '.png']);
print([dataOUTfolder fn '.pdf'],'-dpdf','-bestfit')

%% PLOT TOPO CORR VALUES ALL SUBJECTS
goodsbj = ismember(subjidx, subjidx);
cfg = [];
cfg.layout = layout;
cfg.figure = 'gca';
cfg.commentpos = 'middletop';
cfg.marker = 'on';
cfg.colorbar = 'no';
cfg.style = 'both'; % 'straight';
dd = [];dd.dimord = 'chan_time';dd.label = layout.label;dd.time = 1;
figure(102); set(gcf,'Position',[100 100 900 400]);
p =1;
d = squeeze(mean(gavgr(1).model(goodsbj,:,:),2)) - squeeze(mean(gavgr(3).model(goodsbj,:,:),2)); %
climd = -0.0015;
climu =+0.0015;
% climd = -abs(max([mean(d)]));
% climu = abs(max([mean(d)]));
cfg.zlim = [climd climu];

for i =2:length(models)
    d = squeeze(mean(gavgr(1).model(goodsbj,1:10,:),2)) - squeeze(mean(gavgr(i).model(goodsbj,1:10,:),2)); %
    subplot(2,nmodel-1,p);  dd.avg = mean(d)'; cfg.comment = sprintf(models{i});
    ft_topoplotER(cfg,dd);colormap(flipud(brewermap(64,'RdBu')));

    d = squeeze(mean(gavgr(1).model(goodsbj,11:14,:),2)) -squeeze(mean(gavgr(i).model(goodsbj,11:14,:),2)); %
    subplot(2,nmodel-1,p+nmodel-1);  dd.avg = mean(d)'; cfg.comment = sprintf(models{i});
    ft_topoplotER(cfg,dd);colormap(flipud(brewermap(64,'RdBu')));
    p = p+1;
end
c = colorbar;c.LineWidth = 0.1;
c.Position = [0.92, 0.1, 0.01, 0.08]; % [x, y, width, height]
c.FontSize = 7;
fn = sprintf(['/FIG_GAVG_Rtopo_ALL_DIFF_' ]);
saveas(gcf,[dataOUTfolder fn '.png']);
print([dataOUTfolder fn '.pdf'],'-dpdf','-bestfit')


%% PLOT BARS
k = 1; %ROI
figure(89); set(gcf,'Position',[100 100 1000 800]);
p=1;
for i=2:length(models)
    subplot(2,nmodel-1,p);
    real = mean(r(1).model(goodsbj,1:10, k) ,2) - mean(r(i).model(goodsbj,1:10,k) ,2);
    y = [-0.006 0.008];
    sample_size = numel(real);
    sd = 2*(std(real)/ sqrt(sample_size));
    allmean = [mean(real)];
    bar(allmean);hold on; errorbar(allmean, sd, 'k', 'linestyle', 'none');
    title(sprintf(models{i})); xticklabels({''});ylim(y);
    subplot(2,nmodel-1,p+nmodel-1);
    real = mean(r(1).model(goodsbj,11:14, k) ,2) - mean(r(i).model(goodsbj,11:14, k) ,2);
    sample_size = numel(real);
    sd = 2*(std(real)/ sqrt(sample_size));
    allmean = [mean(real)];
    bar(allmean);hold on; errorbar(allmean, sd, 'k', 'linestyle', 'none');
    title(sprintf(models{i})); xticklabels({''});ylim(y);
    p = p+1;
end
annotation('textbox', [0 0.9 1 0.1], 'String', 'Full Model-each predictor', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14);

fn = sprintf(['/FIG_BAR_ALL_DIFF_' ]);
saveas(gcf,[dataOUTfolder fn '.png']);
print([dataOUTfolder fn '.pdf'],'-dpdf','-bestfit')

%% PLOT WEIGHTs avaeraged across melodies and subjects
time = model(1).mdl{1}{1}.t;
k=1.5;
%plot chn = 34, 53, 17
% goodsbj = ~ismember(subjidx, exclude); %N = 34
avgsubj = squeeze(mean(gavgmodel(1).avg(goodsbj,:,:,:),1));
%avgsubj = squeeze(gavgmodelA(:,:,:,:));

figure(103), set(gcf,'Position',[100 100 800 800]);
Nchan_to_plot = [find(strcmp('Fz', norm_chanlbls)) find(strcmp('CPz', norm_chanlbls)) find(strcmp('Oz', norm_chanlbls))]; %
pred = 1:length(model(1).predname);
p =1;
for n =1:length(Nchan_to_plot)
    chan_to_plot = Nchan_to_plot(n);
    subplot(3,2,p);
    w = zeros(length(pred),length(time));
    w(1:length(pred),:) = mean(avgsubj(pred,:,chan_to_plot),3);xlim([-50 350]);
    if n == 1; limits = [-max(abs(mean(w)))*k max(abs(mean(w)))*k]; end
    clim = limits*k;
    %limits = [min(min(w)) max(max(w))]; ylim(limits);
    imagesc(time,1:size(w,1),w, limits);
    yticks(1:length( model(1).predname)); yticklabels( model(1).predname);
    colorbar('SouthOutside');colormap(flipud(brewermap(64,'PuOr')));
    title(sprintf('Model Full, Chn %s',  norm_chanlbls{chan_to_plot}));
    hold on
    subplot(3,2,p+1);hold on
    for i = 1:size(w,1); plot(time,w(i,:));end
    xlabel('Lags (ms)');ylabel('Model weight');
    if n==3; legend( model(1).predname, 'Location','southeast' );end; ylim([-1.8 1.2]);
    p = p+2;
end
fn = sprintf(['FIG_FULLmodelWEIGHT_' ]);
saveas(gcf,[dataOUTfolder fn '.png']);
print([dataOUTfolder fn '.pdf'],'-dpdf','-bestfit')



%% PLOT CORRELATIONS
songs = 1:10;
figure(90); set(gcf,'Position',[100 100 900 400]);
x = mean(r(1).model(goodsbj,songs) ,2) - mean(r(3).model(goodsbj,songs) ,2);
y = mean(r(1).model(goodsbj,songs) ,2) - mean(r(5).model(goodsbj,songs) ,2);
[rho, pval] = corr(x, y, 'Type', 'Spearman');
subplot(1,2,1); scatter(x, y, 'filled'); hold on;
% Fit a linear model
lm = fitlm(x, y);
% Get the fitted line
fittedX = linspace(min(x), max(x), 100);
fittedY = predict(lm, fittedX');
% Plot the fitted line
plot(fittedX, fittedY, '-r', 'LineWidth', 2);
title(['Spearman : \rho = ', num2str(rho), ', p = ', num2str(pval)]);
xlabel('Delta So');
ylabel('Delta IOI');
grid on;

x = mean(r(1).model(goodsbj,songs) ,2) - mean(r(4).model(goodsbj,songs) ,2);
y = mean(r(1).model(goodsbj,songs) ,2) - mean(r(6).model(goodsbj,songs) ,2);
[rho, pval] = corr(x, y, 'Type', 'Spearman');
subplot(1,2,2); scatter(x, y, 'filled'); hold on;
% Fit a linear model
lm = fitlm(x, y);
% Get the fitted line
fittedX = linspace(min(x), max(x), 100);
fittedY = predict(lm, fittedX');
% Plot the fitted line
plot(fittedX, fittedY, '-r', 'LineWidth', 2);
title(['Spearman : \rho = ', num2str(rho), ', p = ', num2str(pval)]);
xlabel('Delta Sp');
ylabel('Delta ITI');
grid on;

grid on;
sgtitle('REAL condition')
fn = sprintf(['/FIG_Correlation_REAL' ]);
saveas(gcf,[dataOUTfolder fn '.png']);
print([dataOUTfolder fn '.pdf'],'-dpdf','-bestfit')

