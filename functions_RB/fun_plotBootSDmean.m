function [ ] = fun_plotBootSDmean(subjidx,ALLDATA, time_s, chantoplot, xlimits)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
myString = strjoin(cellfun(@string, chantoplot), ',');

for si = 1:length(subjidx)
    [a,ch] = ismember( chantoplot, ALLDATA{si}.label);
    d(si,:,:) = ALLDATA{si}.avg(ch, :);
end


B = 500; bsmean =[]; bsstd=[];stderror=[];
x = d(:,:,:);
x = squeeze(mean(x,2)); %mean over channles
%  Force x to be time*repetitions
x=x';
[mn,sd,bsall] = fBootstrapMean(x,B);
bsmean=mn';
bsstd=sd'; % Bootstrap standard deviation IS the estimate of the sample standard ERROR!
stderror = std(x,0,2)/sqrt(numel(subjidx));


C=[0.7441    0.1754    0.1966];

%%%% Plot boostrap SD and mean across subj for all conditions
% fh =figure(2); clf
% fh.Position = [0, 0, 1200, 500];
subplot(1,4,1);
plot(time_s,bsmean(:),'color',C, 'LineWidth', 1.5); hold on; 
% Plot SE of the bootstrap resampling
b = bsstd'; % STDEV
a = bsmean';
Y = [b+a;  flipud(-b+a)]';
abscissa = time_s';
X = [abscissa; flipud(abscissa)];
h = fill(X,Y,C,'edgecolor','none','facealpha',0.2); hold on;
plot(abscissa,a*0,'-k'); % plot zero line
xlim([xlimits]);
set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
title(['\bf channels: ' myString{1}]); % \bf = bold font
set(gca,'Box','off'); ylabel('\bf Voltage (uV)');


%%% Plot time domain each subject %%%%
subplot(1,4,2);
for si = 1:length(subjidx)
    x = squeeze(d(si,:,:));
    plot(time_s,mean(x,1),'color',C, 'LineWidth', 0.2)  ; hold on
end
plot(abscissa,a*0,'-k');
plot(time_s,bsmean,'color', 'k', 'LineWidth', 1.5);
xlim([xlimits]);
set(gca,'Box','off')
ylabel('\bf Voltage (uV)');
title(['\bf Single subjects']); % \bf = bold fon


x = permute(d,[3 2 1]); % nt_find_outlier_trials requires input structure time * channels * trials)
xb = nt_demean(x, find(time_s<0)); %baseline
 % Reject trials that deviate from the mean by more than twice the average deviation from the mean
good_subj=nt_find_outlier_trials(xb, 2, 0); % creates array with good trials based on threshold
[nout, j]= find(~ismember(1:length(subjidx), good_subj));
d_clean = d(good_subj, :,:);
subplot(1,4,3);
for si = 1:length(good_subj)
    x = squeeze(d_clean(si,:,:));
    plot(time_s,mean(x,1),'color',C, 'LineWidth', 0.2)  ; hold on
end
plot(abscissa,a*0,'-k');
plot(time_s,bsmean,'color', 'k', 'LineWidth', 1.5);
xlim([xlimits]);
set(gca,'Box','off')
ylabel('\bf Voltage (uV)');
title(['\bf removed: ' num2str(j)]); % \bf = bold fon


%%R replot mean  without outliers

B = 500; bsmean =[]; bsstd=[];stderror=[];
x = d_clean(:,:,:);
x = squeeze(mean(x,2)); %mean over channles
%  Force x to be time*repetitions
x=x';
[mn,sd,bsall] = fBootstrapMean(x,B);
bsmean=mn';
bsstd=sd'; % Bootstrap standard deviation IS the estimate of the sample standard ERROR!
stderror = std(x,0,2)/sqrt(numel(subjidx));

subplot(1,4,4);
plot(time_s,bsmean(:),'color',C, 'LineWidth', 1.5); hold on; 
% Plot SE of the bootstrap resampling
b = bsstd'; % STDEV
a = bsmean';
Y = [b+a;  flipud(-b+a)]';
abscissa = time_s';
X = [abscissa; flipud(abscissa)];
h = fill(X,Y,C,'edgecolor','none','facealpha',0.2); hold on;
plot(abscissa,a*0,'-k'); % plot zero line
xlim([xlimits]);
set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
title(['\bf without outliers']); % \bf = bold font
set(gca,'Box','off'); ylabel('\bf Voltage (uV)');