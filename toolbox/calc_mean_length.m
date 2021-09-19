function [km,Sm] = calc_mean_length(S,k)
% This function will compute the energy-weighted inverse length scale
% (note L^-1 = ky/2pi)
% INPUT: S = spectra, k = wavenumbers
% OUTPUT: km = mean inverse length scale, Sm = magnitude of spectra at mean

km  = sum(k.*S)/sum(S); %2pi factor still included
[val,idx] = min(abs(k-km));

if k(idx)-km > 0
    ffrac = (k(idx)-km)/(k(idx)-k(idx-1));
    Sm = ((S(idx-1)-S(idx))*ffrac)+S(idx);
else
    ffrac = (k(idx+1)-km)/(k(idx+1)-k(idx));
    Sm = ((S(idx)-S(idx+1))*ffrac)+S(idx+1);
end
