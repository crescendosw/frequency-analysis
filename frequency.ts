import { PolynomialRegression } from 'ml-regression-polynomial';

// ASDF = Average Squared Difference Function
// Roughly equivalent to an autocorrelation, except here we look for a minimum.
// The minimum corresponds to where a delayed copy of the signal most closely resembles
// the current location of the signal.  This should occur at a distance correlating
// to the fundamental frequency, but should work even when the fundamental frequency
// is missing and we only get the harmonics.
// Input:
//   x : input data vector. Assumed to be longer than kmax
//   kmin,kmax : integer. Together provide the range of offsets to compute.
//       The expected range should be the sampling rate (ie 44kHz) / expected frequency.
//       corresponds to rough_offset
//   out : output data vector with asdf values computed over the range
//         of [kmin,kmax]
export function averageSquaredDifferenceFunction(samples: number[], kmin: number, kmax: number): number[] {
  const {length} = samples;
  if (kmin >= kmax || kmax > length) {
    throw new Error(`Invalid parameter values: sample length = ${length}, kmin = ${kmin}, kmax = ${kmax}`)
  }
  // 1 dimensional array of kmax - kmin + 1 buckets
  // out = zeros(kmax+1-kmin,1);
  const out = Array(kmax - kmin + 1);
  // for k=kmin:kmax
  for (let k = kmin; k <= kmax; k += 1) {  // stepping by 1 (array dereferncing is 1's based)
    // fill sequentially (+1 b/c one based)
    // out(k-kmin+1) = mean((x(1:end-k)-x(k+1:end)).^2);  // array of differences, square each, get the average
    const diffsSquared = Array(length - k - 1);
    for (let i = 0; i < diffsSquared.length; i += 1) {
      const diff = samples[i] - samples[i+k+1];
      diffsSquared[i] = diff * diff;
    }
    out[k - kmin] = mean(diffsSquared);
  }
  return out;
}

// The Newton-Raphson method is a zero-finding algorithm that uses a function and
// its derivative to rapidly converge upon the root.  In general, it is necessary to
// try multiple initial locations (x0) to find all roots.  For our case, it is
// sufficient to try in the vicinity of the expected frequency.
//
// Inputs:
//   f : function with a zero that we wish to locate.
//   df : function that is the derivative of f.
//   x0 : initial estimate of the zero location.
//   tolerance : stop early when |f(x')| < tolerance
//   max_iterations : stop after this many steps.  (Should converge quickly!)
// Outputs:
  //   root : root location such that f(root) ~= 0.

// function root = newton(f, df, x0, tolerance, max_iterations)
export function newtonRaphson(
  f: (x: number) => number,
  df: (x: number) => number,
  x0: number,
  tolerance: number,
  maxIterations: number,
): number {
  //    for i = 1:max_iterations
  for (let i = 0; i < maxIterations; i++) {
    // fx = f(x0);
    const fx = f(x0);
    // dfx = df(x0);
    const dfx = df(x0);

    // if abs(fx) < tolerance
    if (Math.abs(fx) < tolerance) {
      return x0;
    }

    // x0 = x0 - fx / dfx;
    x0 = x0 - fx / dfx;
  }
  return x0;
}

function mean(values: number[]): number {
  return values.reduce((accum, value) => accum + value, 0) / values.length;
}

// Find the error in frequency using the ASDF algorithm.
// ASDF = Average Squared Difference Function
// It is very related to the autocorrelation, except the ASDF is minimized where
// the autocorrelation is maximized.
// Using autocorrelation (or ASDF) is superior to directly using a DFT
// as the autocorelation automatically takes into account harmonics which would
// need to be done manually for the DFT.
// See https://dsp.stackexchange.com/questions/22067/what-is-an-amdf
  // and http://www.musingpaw.com/2012/04/musical-pitch-is-not-just-fft-frequency.html
  // Input is:
  //  x : data samples (in float, about 2048 samples is good choice)
//  freq_expect : expected frequency in Hz
// Output is:
  //  out : frequency error in Hz.
// Note: it may be necessary to use the Log of the error since pitch is on a logrithmic scale.
//
// Example:
//   for C3, f_expect would be 130.81. x would be the most recent set of 2048 samples.
// function out = find_err(x, f_expect)
export function findErrorHz(samples: number[], sampleRate: number, expectedFrequencyHz: number) {
  // Given the expected frequency, compute the rough offset in the ASDF
  // where the minimum occurs.  (The ASDF is equivalent to the autocorrelation,
  // where we would need to find the maximum.)
  //    rough_offset = round(44e3 / f_expect);
  const roughOffset = Math.round(sampleRate / expectedFrequencyHz);

  // Give some bounds to search.  Currently +/- 30 "bins".
  // kmin = rough_offset - 30;
  const kmin = roughOffset - 30;
  // kmax = rough_offset + 30;
  const kmax = roughOffset + 30;

  // Compute the ASDF over that range
  // tmp = asdf(x, kmin, kmax);
  const averageSquaredDifference = averageSquaredDifferenceFunction(samples, kmin, kmax);

  // Fit a 4th order polynomial to the result. We then use the polynomial
  // to search for the minimum.
  // JS library ibrary found in SF article: https://stackoverflow.com/a/73153281/998490
  // Poly = polyfit(kmin:kmax, tmp, 4); // fit 4th order polynomial (5 coefficents)
  const range = Array(kmax - kmin + 1);
  for (let i = kmin; i <= kmax; i += 1) {
    range[i - kmin] = i;
  }
  const degree = 4;
  if (range.length !== averageSquaredDifference.length) {
    throw new Error(
      `range ${range.length} must have same number of values ` +
      `as average squared difference ${averageSquaredDifference.length}`);
  }
  const {coefficients} = new PolynomialRegression(range, averageSquaredDifference, degree);
  coefficients.reverse(); // coefficient order of JS is reversed from Octave

  // der = Poly(1:4).*(4:-1:1);  // take 1st derivative (first 4 coefficents multiplied respecively by [4, 3, 2, 1])
  const derivative = [4 * coefficients[0], 3 * coefficients[1], 2 * coefficients[2], coefficients[3]];
  // der2 = der(1:3).*(3:-1:1);  // take 2nd derivative (first 3 coefficents multiplied respecively by [3, 2, 1])
  const secondDerivative = [3 * derivative[0], 2 * derivative[1], derivative[2]];

  // Minimization: set the derivative to 0 and solve for x.
  // We can use the Netwon-Rhapson method for this since
  // we can easily compute the derivative and it's (the 2nd) derivative.
  // create some lambda functions to pass in to the Newton-Rhapson.
  // polyval = given an x, evaluate the coefficents (see Brian's notes on not losing precision)
  //    f = @(x)polyval(Poly,x); df = @(x)polyval(der,x); df2 = @(x)polyval(der2,x);
  const df = (x: number) => ((derivative[0] * x + derivative[1]) * x + derivative[2]) * x + derivative[3];
  const df2 = (x: number) => (secondDerivative[0] * x + secondDerivative[1]) * x + secondDerivative[2];

  // Call the NR (Newton-Rhapson) method to get the minimizing offset
  // off = newton(df,df2, rough_offset, 1e-10, 30); // 30 unrelated to 30 above (max iterations)
  const off = newtonRaphson(df, df2, roughOffset, 0.0000000001, 30);

  // There are other ways to obtain roots.  This one is fairly efficient in part
  // because it is only looking in the region of interest to us.  Because the
  // 3rd order polynomial (der) is pretty smooth, this algorithm also converges
  // quickly.

  // Convert from offset to frequency, assuming 44kHz sampling
  // f_est = 44e3 / off;
  const frequency = sampleRate / off;

  // Compute the error relative to the expected frequency
  // out = f_est - f_expect;
  return frequency - expectedFrequencyHz;
}
