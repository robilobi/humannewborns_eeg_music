function endtriggers = RB_getsonglength(wav_path, stim_names)
% return song duration in sec, the order of the songs is the one given as
% input
endtriggers = [];
for f = 1:length(stim_names)
    disp(stim_names{f}); % show which file is being loaded
    % Load the sound file
    [y,stimFs] = audioread([wav_path stim_names{f} '.wav']);
    duration = length(y) / stimFs;
    endtriggers = [endtriggers round(duration, 3)];
end
