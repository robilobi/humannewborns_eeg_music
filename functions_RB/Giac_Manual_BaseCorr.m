function [ data_out ] = Giac_Manual_BaseCorr( data, interval )
%%
% Script doing baseline correction 'manually'. Data is the input
% structure. Interval express time range of the baseline (in msec!). The
% advantage over using a canonical fieldtrip structure is that data keep
% the original format and is not changed into a rpt_chan_time format.
%% Giacomo Novembre

    data_out = data;

    for t = 1:length(data.trial)
        for c = 1:length(data.label)

            tr_tmp      = data.trial{t}(c,:);
            tr_tmp_time = data.time{t};
            
                %str = find(tr_tmp_time==interval(1)./1000);
                %nde = find(tr_tmp_time==interval(2)./1000);

                [~,closest_start] = min(abs(tr_tmp_time-(interval(1)./1000))); % find values that are closest to the desired ones
                [~,closest_end]   = min(abs(tr_tmp_time-(interval(2)./1000)));
                
                    baseline = mean(tr_tmp(closest_start:closest_end));
                    tr_tmp_bc = tr_tmp - baseline;
                    
                        data_out.trial{t}(c,:) = tr_tmp_bc;
                    
        end
    end

end

