% ADD IdYOM preditocrs to stim Matrix
clear all
%%% PATHS
U = 'C:\Users\robianco\OneDrive - Fondazione Istituto Italiano Tecnologia';
addpath('.\functions_RB')
addpath([U '\MATools\mTRF-Toolbox-master\mtrf']); % add the mTRF toolbox functions

stim_path = ('..\STIMULI\');
midi_path = ('..\STIMULI\midiall\');
dataINFolder = ('..\RESULTS\');

%%% VIEWPOINTS: CPITCH IOI-RATIO
% idyomfl(1).name = '89-onset-ioi-ratio-nil-nil-melody-nil-10-both+-nil-t-nil-c-nil-t-t-x-3.dat';
% idyomfl(2).name = '89-cpitch-cpitch-nil-nil-melody-nil-10-both+-nil-t-nil-c-nil-t-t-x-3.dat';
% outname = 'FINAL_StimMatrix_both+';

idyomfl(1).name = '89-onset-ioi-ratio-1_2_3-nil-melody-nil-10-both-nil-t-nil-c-nil-t-t-x-3.dat';
idyomfl(2).name = '89-cpitch-cpitch-1_2_3-nil-melody-nil-10-both-nil-t-nil-c-nil-t-t-x-3.dat';
outname = 'FINAL_StimMatrix_both';


namecut = 2; % = 2 for 89; = 4 for idyom 93
% Load the amplitude envelopes spectral flux
predictor = 'StimAcoustics'; %
allpred = load([stim_path predictor]);

Fs = 100; % must be the same of the eeg and the acoustic descriptors
sound_delay = 5; % delay of the actual stimulus start relative to the start of the EEG trial
% (needed to adjust the onset times of the IDyOM model values)
% (20-4-2022) The stimulus onset  is 5 s after relative to the
% beep at the beginning of the .wav file
% There is also a beep at the end of the .wav file. The time from
% stimulus end to the beep *offset* is to be 5 s as well

midi_names = {'audio01.mid','audio02.mid','audio03.mid','audio04.mid','audio05.mid', ...
    'audio06.mid','audio07.mid','audio08.mid', 'audio09.mid','audio10.mid',...
    'shf01.mid','shf05.mid','shf08.mid','shf10.mid'};

for s = 1:length(midi_names)  
    filename = midi_names{s}(1:end-4);
    disp(filename); % show which file is being loaded
    %find match with acoustic
    pred_idx = find(strcmp(allpred.stim_names,[filename '.wav']));
    cutpred = allpred.spflux{pred_idx}; %take max length of predictor

    mid = readmidi([midi_path '/' filename '.mid']);
    noteinfo = midiInfo(mid,0);
    pitch = noteinfo(:,3);
    time_onsets = noteinfo(:,5);

    % ONSET
    zero_pad_end=3*sound_delay;
    dur = sound_delay+max(time_onsets)+zero_pad_end; % add 1 second after the last no
    feat=delta_vec(time_onsets,ones(length(time_onsets),1),sound_delay,dur,Fs);
    allpred.onset{s} = feat(1:length(cutpred));

    % ITI
    iti = [99 ; diff(pitch)];
    iti = abs(iti);
    feat = delta_vec(time_onsets,iti,sound_delay,dur,Fs);
    allpred.iti{s} = feat(1:length(cutpred));

    % IOI
    ioi = [99 ; diff(time_onsets)];
    feat = delta_vec(time_onsets,ioi,sound_delay,dur,Fs);
    allpred.ioi{s} = feat(1:length(cutpred));

    %SO ioi-ratio
    headers = {'melody.name', 'onset.information.content', 'onset.entropy'};
    Sonset = Nate_get_idyom_vars([stim_path idyomfl(1).name],headers,[0 1 1 ]);
    Sonset.melody_name = cellfun(@(x) x(namecut:end-1), Sonset.melody_name, 'UniformOutput', false); %make name as  midi_names
    stim_idx = cellfun(@(x) strcmp(x,filename), Sonset.melody_name);
    So = Sonset.onset_information_content(stim_idx);
    feat = delta_vec(time_onsets,So,sound_delay,dur,Fs);
    allpred.Soioiratio{s} = feat(1:length(cutpred));
    %EO
    stim_idx = cellfun(@(x) strcmp(x,filename), Sonset.melody_name);
    Eo = Sonset.onset_entropy(stim_idx);
    feat = delta_vec(time_onsets,Eo,sound_delay,dur,Fs);
    allpred.Eoioiratio{s} = feat(1:length(cutpred));


    %SP pitch
    headers = {'melody.name', 'cpitch.information.content', 'cpitch.entropy'};
    Spitch = Nate_get_idyom_vars([stim_path idyomfl(2).name],headers,[0 1 1 ]);
    Spitch.melody_name = cellfun(@(x) x(namecut:end-1), Spitch.melody_name, 'UniformOutput', false); %make name as  midi_names
    stim_idx = cellfun(@(x) strcmp(x,filename), Spitch.melody_name);
    Sp = Spitch.cpitch_information_content(stim_idx);
    feat = delta_vec(time_onsets,Sp,sound_delay,dur,Fs);
    allpred.Sppitch{s} = feat(1:length(cutpred));
    %EP
    stim_idx = cellfun(@(x) strcmp(x,filename), Spitch.melody_name);
    Ep = Spitch.cpitch_entropy(stim_idx);
    feat = delta_vec(time_onsets,Ep,sound_delay,dur,Fs);
    allpred.Eppitch{s} = feat(1:length(cutpred));

    allpred.env{s} = allpred.env{s}(1:length(cutpred));

    stim_all{s} = [allpred.spflux{s}, allpred.onset{s}, allpred.iti{s}, allpred.ioi{s}, ...
                    allpred.Soioiratio{s}, allpred.Eoioiratio{s}, ...
                    allpred.Sppitch{s}, allpred.Eppitch{s}];

end

pred_names = {'Spflux', 'onset', 'ITI', 'IOI', ...
               'Soioiratio', 'Eoioiratio', ...
               'Sppitch','Eppitch' };

save([stim_path outname], 'stim_all', 'pred_names');

%% [all_id pitch onsets alliti allioi -So Eo - Sp Ep;
load([stim_path outname])
% it assumes melodies are in the order of midi_names
allinfo = [];
for i = 1:length(stim_all)
    filename = midi_names{i}(1:end-4);
    disp(filename); % show which file is being loaded
    mid = readmidi([midi_path '/' filename '.mid']);
    noteinfo = midiInfo(mid,0);
    pitch = noteinfo(:,3);
    time_onsets = noteinfo(:,5);
    thismel = stim_all{i};
    idx = find(thismel(:,2)>0);
    melid = repmat(i, length(idx),1);
    melinfo     = [melid  pitch time_onsets thismel(idx,3) thismel(idx,4) ...
                    thismel(idx,5) thismel(idx,6) ...
                    thismel(idx,7) thismel(idx,8)];
    allinfo =[allinfo; melinfo];
end
writematrix(allinfo, [stim_path 'allinfo_stimuli_' outname '.csv']);

%% plot stim predictors
i = 4;
thismel = stim_all{i};
use_idx = (5.03*100+1):(5.03*100+300);
thismel =thismel(use_idx,:);
thismel =thismel./(ones(length(use_idx),1)*rms(thismel));
count = 1;
for j = [2 3 4 5 8]
    subplot(5,1,count);
    plot(thismel(:,j))
    ylim([0 10]);
    count = count +1;
end
print([stim_path 'predictors.pdf'],'-dpdf','-bestfit')

%% Functions %%
function dlt = delta_vec(t,vals,sound_delay,dur,Fs)
% Create a delta function with non-zero values at times t
dlt = zeros(ceil(dur*Fs),1);
idx = round((t+sound_delay)*Fs)+1;
dlt(idx) = vals;
end