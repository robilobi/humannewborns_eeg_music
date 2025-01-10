function [ data_dss ] = fun_DSSonset( data, comptokeep, latencydss, referencezero)
% Roberta Bianco 2020
% DSS function as demontrasted by http://audition.ens.fr/adc/NoiseTools/
% DSS finds the most repeatable component across all trials based on a
% specific time window as weighted bias matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v
%referencezero = 0;  % from where to compute the baseline for dss

%% 1. prepare data
datamat = nt_trial2mat(data.trial); % creates structure time * channel * trials
time = data.time{1};
data_dss = data;

[val, time_index(1)] = min(abs(time-latencydss(1))); % returns cell corresponding to t = 0
[val, time_index(2)] = min(abs(time-latencydss(2))); % returns cell corresponding to t = 0.4
time2 = time(time_index(1): time_index(2)); % array of time steps between 0 - 0.4
x = datamat(time_index(1):time_index(2),:,:); % new data structure from datamat with above time interval
% baseline/demean
%x = nt_demean(x); % use nt_demean() instead if there's a strong trend
x = nt_demean2(x, find(time2<referencezero)); % x is a 3D matrix time x channels x trials

%% 2. do DSS
[dssweights,pwr0,pwr1] = nt_dss1(x); % todss = 'denoising matrix', i.e. nchanns x nchanns matrix of weights for each component
comp = nt_mmat(datamat,dssweights); % matrix multiplication of data structure and denoising matrix. time * components * trial.
f1 = figure(1); clf; plot(pwr1./pwr0,'.-'); ylabel('score'); xlabel('component');

f2 = figure(2); clf;
for iComp=1:6
    subplot(2,5,iComp);
    nt_bsplot(comp(:,iComp,:),[],[],time);
    title(iComp);
    xlim([time(1) time(end)]); xlabel('s');
end

%% 3. Project back into sensor space the components that you want
% (do this after you have run the above dss on all subjects to extract & plot the compone t timeseries, and decided how
% many components to keep!!)
KEEP = 1:comptokeep; % all of the components to keep, i.e. if you want first 3, enter KEEP = 1:3
%ccov = nt_xcov(comp,datamat); % c is cross-covariance between component and data across all time steps. Components * channel.
ccov = nt_xcov(comp, datamat)/(size(comp,1)*size(comp,3));
datamat_dssed = nt_mmat(comp(:,KEEP,:),ccov(KEEP,:)); % matrix multiplication of kept components of ccv and comp matrix, creating new time * channel * trial matrix
data_dss.trial = nt_mat2trial(datamat_dssed); % replaces data.trial with dssed data

end
