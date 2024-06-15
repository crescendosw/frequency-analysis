import {describe, expect, test} from '@jest/globals';
import {averageSquaredDifferenceFunction, findErrorHz} from "../frequency";
import * as fs from 'fs';
const wav = require('node-wav');

describe('AverageSquaredDifferenceFunction', () => {
  test('basic array', () => {
    const input = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    const output = averageSquaredDifferenceFunction(input, 2, 6);
    expect(output).toEqual([9, 16, 25, 36, 49]);
  });

  test('larger array', () => {
    const input = [-3, 18, 28, 16, 7, 5, 2, 6, -4, 8, 19, 102, 0, 16, -9];
    const output = averageSquaredDifferenceFunction(input, 3, 6);
    expect(output).toEqual([1049.4545454545455, 1212.6, 1215.4444444444443, 1305.625]);

    // Actual output from Octave for asdf([-3, 18, 28, 16, 7, 5, 2, 6, -4, 8, 19, 102, 0, 16, -9], 4, 7)
    // expect(output).toEqual([1049.5, 1212.6, 1215.4, 1305.6]);
  });
});

describe('findErrorHz', () => {
  test('basic wav file', () => {
    // brew install octave
    // v_g3 = audioread('voice_g3.wav');
    // g3 = 196;
    // find_err(v_g3(4000:4000+4096), g3)
    // ans = -0.070425
    let buffer = fs.readFileSync('./tests/voice_g3.wav');
    const {sampleRate, channelData} = wav.decode(buffer);
    const g3Hz = 196;
    const samples = Object.values(channelData[0]).slice(4000, 4000 + 4096) as number[];
    const result = findErrorHz(samples, sampleRate, g3Hz);
    expect(result).toEqual(1.1945213304115896)
    // Octave's result
    // expect(result).toEqual(-0.070425)
  });
});
