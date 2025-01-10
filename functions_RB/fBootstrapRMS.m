function [a,b,all]=fBootstrapRMS(x,N)
%[rms,sdev,all] = fBoostrapRMS(x,N) - calculate rms, estimate sd using bootstrap
%
%  a: rms over second dimension of x
%  b: standard deviation of rms calculated by bootstrap
%  all: matrix of all boostrap iterations
%
%  x: matrix of observations (time X repetitions)
%  B: number of bootstrap trials [default: 500]
%
% Nicolas 15/07/2015

if nargin <2; N=500; end

if ~ismatrix(x); error('Input matrix must be at most 2D'); end

[nsamples,nrep]=size(x);
all=zeros(nsamples,N);
for k=1:N
    % Choose random combination of indices (with replacement)
    idx=randi(nrep,1,nrep);
    all(:,k)=sqrt(mean(x(:,idx).^2,2));
end

a=rms(x,2);
b = nanstd(x',1)/sqrt(nrep);
%b = nanstd(x',1);
%b=b';
%b=rms(all-repmat(a,1,N));
%b=sqrt(mean((all-repmat(a,1,N)).^2,2));

plot(x)
% h=[];
% h = [h;plot(time,a,'LineWidth',1.2)]; hold on;
% Y = [b+a;flipud(-b+a)]';
% X = [time(:); flipud(time(:))];
% hp = fill(X,Y, 'b','edgecolor','none','facealpha',0.2); hold on;
% set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
% ylim([0 1.2e-12]);


%se = nanstd(all',1);
%sdev=sqrt(mean((all-repmat(rmsout,1,B)).^2,2));