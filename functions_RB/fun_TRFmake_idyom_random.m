function [idyom_shu] = fun_TRFmake_idyom_random(idyom_v, var_to_randomise)
%(spflux, ONSET, 'ITI','IOI','Sp','So', 'Ep', 'Eo')

idyom_shu = idyom_v;
idx = find(idyom_v(:,2) > 0); % where onset is 1

for i = 1:length(var_to_randomise)
    col =  var_to_randomise(i);
    rng(i)
    idxshu = idx(randperm(length(idx)));
    idyom_shu(idx, col) = idyom_v(idxshu, col);
end
end