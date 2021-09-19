function [emean] = calc_energy_weighted_mean(spec,freq,freqrange)
% calc energy weighted mean

[val, idmin] = nanmin(abs(freq-freqrange(1)));
[val, idmax] = nanmin(abs(freq-freqrange(2)));
specint = trapz(freq(idmin:idmax),spec(idmin:idmax));
emean = specint/(freq(idmax)-freq(idmin));