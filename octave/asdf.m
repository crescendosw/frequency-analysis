% ASDF = Average Squared Difference Function
% Roughly equivalent to an autocorrelation, except here we look for a minimum.
% The minimum corresponds to where a delayed copy of the signal most closely resembles
% the current location of the signal.  This should occur at a distance correlating
% to the fundamental frequency, but should work even when the fundamental frequency 
% is missing and we only get the harmonics.
% Input:
%   x : input data vector. Assumed to be longer than kmax
%   kmin,kmax : integer. Together provide the range of offsets to compute.
%       The expected range should be the sampling rate (ie 44kHz) / expected frequency.
%   out : output data vector with asdf values computed over the range
%         of [kmin,kmax]
function out = asdf(x, kmin=1, kmax=200)
   out = zeros(kmax+1-kmin,1);
   for k=kmin:kmax
      out(k-kmin+1) = mean((x(1:end-k)-x(k+1:end)).^2);
end
