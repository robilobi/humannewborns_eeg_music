function [pos, neg]=fClusterbasedPermTest(conds1, conds2, N, time, clusteralpha)

            % perform the statistical test using randomization and a clustering approach
            cfg = [];
            %cfg.avgoverchan      = 1;
            cfg.statistic        = 'ft_statfun_depsamplesT';
            cfg.latency          = [time(1) time(end)];
            cfg.numrandomization = 1000;
            cfg.clusterthreshold = 'nonparametric_common';
            cfg.correctm         = 'cluster';  %'cluster';
            cfg.method           = 'ft_statistics_montecarlo';  % use the Monte Carlo method to calculate probabilities
            cfg.tail             = 0;  % Set to 0 for two-tailed, 1 or -1 for one-tailed tests
            cfg.clusteralpha     = clusteralpha;
            cfg.alpha            = 0.05;

            % Sample dependent samples design for two conditions and multiple trials/subjects
            num_subjects = N;  % number of trials/subjects per condition
            cfg.design = [1:num_subjects, 1:num_subjects;  % subject indices
                ones(1, num_subjects), 2*ones(1, num_subjects)];  % condition labels
            cfg.uvar = 1;  % unit of observation (subjects)
            cfg.ivar = 2;  % independent variable (condition)
            cfg.dimord = 'time_subj';
            cfg.dim=[1,numel(time)];
          
            stat = ft_statistics_montecarlo(cfg, [conds1 conds2],cfg.design);

            % Find indices of significant clusters
            pos=[]; neg=[];
            if isfield(stat,'posclusters')
                if ~isempty(stat.posclusters)
                    pos_cluster_pvals = [stat.posclusters(:).prob];
                    pos_signif_clust = find(pos_cluster_pvals < cfg.alpha);
                    poss = squeeze(ismember(stat.posclusterslabelmat, pos_signif_clust));
                    pos = [find(diff([0; poss; 0])==1) find(diff([0; poss; 0])==-1)];
                   if  ~isempty(pos); if pos(end) > size(time); pos(end) = pos(end)-1; end; end

                end
            end
            if isfield(stat, 'negclusters')
                if ~isempty(stat.negclusters)
                    neg_cluster_pvals = [stat.negclusters(:).prob];
                    neg_signif_clust = find(neg_cluster_pvals <cfg.alpha);
                    negs = squeeze(ismember(stat.negclusterslabelmat, neg_signif_clust));
                    neg = [find(diff([0; negs; 0])==1) find(diff([0; negs; 0])==-1)];
                    if  ~isempty(neg)
                    if neg(end) > size(time); neg(end) = neg(end)-1; end
                    end

                end
            end
            
