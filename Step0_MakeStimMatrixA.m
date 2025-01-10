% Compute the envelopes / Env Derivative / Spectral flux
%%% COMPUTE ENVELOPE AND DERIVATIVE AS DONE IN E. LALOR'S LAB
% load functions for envelope extraction
U = 'C:\Users\robianco\OneDrive - Fondazione Istituto Italiano Tecnologia\';
addpath(genpath([U '\MATools/Envelope Extraction/']));
%addpath([U '\MATools/NoiseTools']);

Fs = 100; % final sampling rate of the envelope !!! 100 it's default srate of spectral flux output
stim_paths = ([U '\BACHNB\STIMULI\test_wav_final_beep\']);
dataOUTfolder = [U '\BACHNB\STIMULI\']; mkdir(dataOUTfolder);

stim_names = {};
env = {};
denv = {};
spflux = {};
disp(stim_paths);

% get the files contained in those paths
fls = dir([stim_paths]);
for ii = 1:length(fls)
    % if it's a .wav file
    if length(fls(ii).name)>4
        if strcmp(fls(ii).name(end-3:end),'.wav')
            disp(fls(ii).name); % show which file is being loaded
            % Load the sound file
            [y,stimFs] = audioread([stim_paths '/' fls(ii).name]);
            % Calculate the amplitude envelope
            eg = gammatoneenv(mean(y,2),stimFs)';
            e = abs(hilbert(mean(y,2)));
            % Downsample the envelope
            eg = nt_resample(eg,Fs,stimFs);
            e = nt_resample(e,Fs,stimFs);

            plot(zscore(eg)); hold on;  plot(zscore(e));
            % Calculate the derivative of the envelope (half-wave
            % rectified)
            de = [0; diff(e)];
            de(de<0) = 0;

            hold on; plot(zscore(de));
            % store the envelope and the derivative
            env = [env; {e}];
            denv = [denv; {de}];

            [s,cf,t] = melSpectrogram(y,stimFs);
            %  imagesc(1:length(s),cf, s);
            s = squeeze(s(:, :, 1));             % take only one channel
            % f = spectralFlux(s,cf);    %% this should be done on the same sample rate on the derivative of env
            currentfs = round(size(s,2)/t(end)); % compute the current sampling rate (it should be higher that the desired final one)
            sresampled = nt_resample(s',Fs, currentfs);
            % imagesc(sresampled');
            f = [spectralFlux(sresampled',cf); 0]; %this should be done on the same samp rate on the derivative of env
            hold on; plot(zscore(f));

            spflux = [spflux; {f(:,1)}];
            % store the stimulus name
            stim_names = [stim_names; {fls(ii).name}];


            %             % Calculate the spectrogram
            [sgam,GCresp] = gammatonespec(mean(y,2),stimFs,32);
            sgam = sgam';
            % Downsample the spectrogram
            sgam = nt_resample(sgam,Fs,stimFs);
            sgam(sgam<0)=0;
            % Calculate the derivative : spectral flux
            sflux = [zeros(1,32) ; diff(sgam)];
            sflux(sflux<0) = 0; %(edited)
            sp = mean(sflux, 2);
        end
    end
end

% Save the stimulus envelopes
save([dataOUTfolder 'StimAcoustics'],'stim_names','env','denv', 'spflux');


