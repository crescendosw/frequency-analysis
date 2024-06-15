% Find the error in frequency using the ASDF algorithm.
% ASDF = Average Squared Difference Function
% It is very related to the autocorrelation, except the ASDF is minimized where
% the autocorrelation is maximized.
% Using autocorrelation (or ASDF) is superior to directly using a DFT
% as the autocorelation automatically takes in to account harmonics which would
% need to be done manually for the DFT.
% See https://dsp.stackexchange.com/questions/22067/what-is-an-amdf
% and http://www.musingpaw.com/2012/04/musical-pitch-is-not-just-fft-frequency.html
% Input is:
%  x : data samples (in float, about 2048 samples is good choice)
%  freq_expect : expected frequency in Hz
% Output is:
%  out : frequency error in Hz.
% Note: it may be necessary to use the Log of the error since pitch is on a logrithmic scale.
%
% Example: 
%   for C3, f_expect would be 130.81. x would be the most recent set of 2048 samples.
function out = find_err(x, f_expect)
   % Given the expected frequency, compute the rough offset in the ASDF
   % where the minimum occurs.  (The ASDF is equivalent to the autocorrelation,
   % where we would need to find the maximum.)
   rough_offset = round(44e3 / f_expect);
   % Give some bounds to search.  Currently +/- 30 "bins".
   kmin = rough_offset - 30;
   kmax = rough_offset + 30;

   % Compute the ASDF over that range
   tmp = asdf(x, kmin, kmax);

   % Fit a 4th order polynomial to the result. We then use the polynomial
   % to search for the minimum.
   Poly = polyfit(kmin:kmax, tmp, 4); % fit 4th order polynomial
   der = Poly(1:4).*(4:-1:1);  % take 1st derivative
   der2 = der(1:3).*(3:-1:1);  % take 2nd derivative

   % Minimization: set the derivative to 0 and solve for x.
   % We can use the Netwon-Rhapson method for this since
   % we can easily compute the derivative and it's (the 2nd) derivative.
   % create some lambda functions to pass in to the Newton-Rhapson.
   f = @(x)polyval(Poly,x); df = @(x)polyval(der,x); df2 = @(x)polyval(der2,x);
   % Call the NR method to get the minimizing offset
   off = newton(df,df2, rough_offset, 1e-10, 30);
   % There are other ways to obtain roots.  This one is fairly efficient in part
   % because it is only looking in the region of interest to us.  Because the
   % 3rd order polynomial (der) is pretty smooth, this algorithm also converges
   % quickly.

   % Convert from offset to frequency, assuming 44kHz sampling
   f_est = 44e3 / off;
   % Compute the error relative to the expected frequency
   out = f_est - f_expect;
   return
end
