function [ data_out ] = Giac_2d_3d_ft_Converter( data, task )
%
% Turn Fieldtrip data structure from 2d to 3d or vice versa. Data is the
% data structure (input), task specifies the type of conversion to be done
% (2d_to_3d or 3d_to_2d). The function assumes that all trials have the
% same duration (if they don't, all trials are concatenated into 1 single
% trial).
%
%% Giacomo Novembre

data_out = data;

switch task
    case '2d_to_3d'
       
        trials_dur = [];
        for kk=1:length(data.trial)
        trials_dur = [trials_dur length(data.trial{kk})]; % compute duration of all trials
        end
        
        if length(unique(trials_dur)) == 1 % if all trials have equal duration, keep trials
        
            new_tr_str = nan([length(data.trial) size(data.trial{1,1})]); 
                for tr = 1:length(data.trial)  
                    new_tr_str(tr,:,:) = data.trial{1,tr};   
                end
            time_epoch = data.time{1,tr};
            
        elseif length(unique(trials_dur)) > 1 % if all trials have NOT equal duration, concatenate trials
            
            new_tr_str = nan(1, length(data.label), sum(trials_dur));
            new_tr_str(1,:,:) = cat(2,data.trial{:});
            time_epoch = 1:sum(trials_dur);
            
        end
        
        data_out.trial  = new_tr_str;
        data_out.time   = time_epoch;
        data_out.dimord = 'rpt_chan_time';
            
    case '3d_to_2d'

        if strcmp('rpt_chan_time',data.dimord)==1 % preliminary check
            tr_nr = size(data.trial,1);
        else
            display('GIAC: dimension error!');
            stop
        end
        
        d = struct;
        d.trial = cell(1,tr_nr);
        d.time  = cell(1,tr_nr);

        for tr = 1:tr_nr
            d.trial{1,tr} = squeeze(data.trial(tr,:,:));
            d.time{1,tr}  = data.time; 
        end
   
        data_out.trial = d.trial;
        data_out.time  = d.time;
        
        data_out = rmfield(data_out,'dimord');
        
end

end

