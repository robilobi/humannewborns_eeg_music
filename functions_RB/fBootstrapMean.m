function [m,sdev,all]=fBootstrapMean(x,B)
%% [m,sd,all]=bsmean(x,B) - calculate mean, estimate sd using bootstrap
% Output:
%  mn: mean of x over second dimension
%  sd: standard deviation from mn of bootstrap trials
%  all: matrix of all bootstrap trials
% Input: 
%  x: matrix of observations (time X repetitions)
%  B: number of bootstrap trials
%
% Nicolas 15/07/2015

if nargin <2; B = 1000; end 

if ~ismatrix(x)
    x=squeeze(x);
    if ~ismatrix(x); error('data must be at most 2D'); end
end

[m,n]=size(x);
all=zeros(m,B);
for k=1:B
    idx=ceil(n*rand(1,n));
    all(:,k)=nanmean(x(:,idx),2);
end

m=nanmean(x,2);
sdev=sqrt(nanmean((all-repmat(m,1,B)).^2,2));