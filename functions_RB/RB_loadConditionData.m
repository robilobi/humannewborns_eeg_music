function data_out = RB_loadConditionData(filename,prestimulus, durations, triggers, trg)
%% load and cut condition data
% insert name of the file: filename, e.g., 'XXX.vhdr'
% define how much to cut before trigger (in sec) after trigger (in sec, customise after trigger)
% eventvalue, e.g., 'S  1' in cell
% eventvalue, e.g., 1 in numbers
% RB (nov 22)

% load the data
cfg                     = [];
cfg.dataset             = filename;
data_proc               = ft_preprocessing(cfg);

% extract event structure
event = ft_read_event(filename);
event = struct2table(event);

cfgtr =[];
count = 1;
for i = 1:height(event)
    idx = find(strcmp(event.value(i), triggers));
    if idx > 0
        cfgtr.trl(count,1) = event.sample(i) - prestimulus*data_proc.fsample;                  % onset
        cfgtr.trl(count,2) = event.sample(i) + floor(durations(idx)*data_proc.fsample);        % offset
        cfgtr.trl(count,3) = -prestimulus*data_proc.fsample;                                                                % trigger time
        cfgtr.trl(count,4) = trg(idx);
        count = count +1;
    end
end

data_out = ft_redefinetrial(cfgtr, data_proc);

end