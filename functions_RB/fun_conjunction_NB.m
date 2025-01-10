function [values, order_chan] = fun_conjunction_NB(stats, filename)
load('Layout_Newborn_EEG.mat');
order_chan = lay.label;
h=figure;clf
h.Position = [100 100 800 700];
time = stats.time;
m1_label = stats.label;
[a, ord] = ismember(order_chan, m1_label);
m1 = stats.mask(ord,:);
mconj = m1;
mconj(mconj~=2) = 0; %set non common values to zero
mconj(mconj==2) = 1; %set non common values to zero

for s = 1:2
    limits = [-6 6];
    stat = stats{1,s};
    if ~isfield(stat, 'posclusterslabelmat')
        stat.posclusterslabelmat = zeros(size(order_chan,1),length(time));end
    if ~isfield(stat, 'negclusterslabelmat')
        stat.negclusterslabelmat = zeros(size(order_chan,1),length(time));end
    M = stat.mask(order(:,s),:).*(stat.posclusterslabelmat(order(:,s),:)-stat.negclusterslabelmat(order(:,s),:));
    subplot(1,3,s);
    imagesc(time,1:length(order_chan),M, limits);hold on;
    yticks(1:size(order_chan,1));
    yticklabels(order_chan);
    title(['Mk ' num2str(s)])
    colormap(flipud(brewermap(64,'RdBu')));
end

subplot(1,3,s+1);
imagesc(time,1:length(order_chan),M.*mconj, limits);
yticks(1:size(order_chan,1));
yticklabels(order_chan);title('Conjunction');
saveas(gcf, [filename '.png']);
print([filename '.pdf'],'-dpdf','-bestfit')


values  = M.*mconj;

end