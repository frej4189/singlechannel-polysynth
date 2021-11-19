//TODO: PHASE SHIFTING

#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>
#include <complex>

using namespace std;

const double C0=16.35;
const double DB0=17.32;
const double D0=18.35;
const double EB0=19.45;
const double E0=20.6;
const double F0=21.83;
const double GB0=23.12;
const double G0=24.5;
const double AB0=25.96;
const double A0=27.5;
const double BB0=29.14;
const double B0=30.87;
const double C1=32.7;
const double DB1=34.65;
const double D1=36.71;
const double EB1=38.89;
const double E1=41.2;
const double F1=43.65;
const double GB1=46.25;
const double G1=49;
const double AB1=51.91;
const double A1=55;
const double BB1=58.27;
const double B1=61.74;
const double C2=65.41;
const double DB2=69.3;
const double D2=73.42;
const double EB2=77.78;
const double E2=82.41;
const double F2=87.31;
const double GB2=92.5;
const double G2=98;
const double AB2=103.83;
const double A2=110;
const double BB2=116.54;
const double B2=123.47;
const double C3=130.81;
const double DB3=138.59;
const double D3=146.83;
const double EB3=155.56;
const double E3=164.81;
const double F3=174.61;
const double GB3=185;
const double G3=196;
const double AB3=207.65;
const double A3=220;
const double BB3=233.08;
const double B3=246.94;
const double C4=261.63;
const double DB4=277.18;
const double D4=293.66;
const double EB4=311.13;
const double E4=329.63;
const double F4=349.23;
const double GB4=369.99;
const double G4=392;
const double AB4=415.3;
const double A4=440;
const double BB4=466.16;
const double B4=493.88;
const double C5=523.25;
const double DB5=554.37;
const double D5=587.33;
const double EB5=622.25;
const double E5=659.25;
const double F5=698.46;
const double GB5=739.99;
const double G5=783.99;
const double AB5=830.61;
const double A5=880;
const double BB5=932.33;
const double B5=987.77;
const double C6=1046.5;
const double DB6=1108.73;
const double D6=1174.66;
const double EB6=1244.51;
const double E6=1318.51;
const double F6=1396.91;
const double GB6=1479.98;
const double G6=1567.98;
const double AB6=1661.22;
const double A6=1760;
const double BB6=1864.66;
const double B6=1975.53;
const double C7=2093;
const double DB7=2217.46;
const double D7=2349.32;
const double EB7=2489.02;
const double E7=2637.02;
const double F7=2793.83;
const double GB7=2959.96;
const double G7=3135.96;
const double AB7=3322.44;
const double A7=3520;
const double BB7=3729.31;
const double B7=3951.07;
const double C8=4186.01;
const double DB8=4434.92;
const double D8=4698.63;
const double EB8=4978.03;
const double E8=5274.04;
const double F8=5587.65;
const double GB8=5919.91;
const double AB8=6644.88;
const double A8=7040;
const double BB8=7458.62;
const double B8=7902.13;

double sine(double time, double freq, double phase) {
    return sin(freq * time + phase);
}

double square(double time, double freq, double phase) {
    return round(sine(time, freq, phase));
}

double sawtooth(double time, double freq) {
    double period = 1.0/freq;
    double result = (fmod(time, period) / period);
    //cout << time << " " << period << " " << result << endl;
    return result;
}

double triangle(double time, double freq, double phase, double last) {
    double sinval = square(time, freq, phase);
    return sinval + last;
}


enum WaveType {SINE, SQUARE, SAWTOOTH, TRIANGLE};

class Oscillator {
    public:
    WaveType waveType;
    int sampleSize;
    double time = 0.0;
    double counter = 0;
    double rps = 0;
    long MAX_AMPLITUDE;
    double val = 0;
    double phase = 0.0;
    vector<long> stream;

    Oscillator(WaveType wave, int sample, int depth) {
        waveType = wave;
        sampleSize = sample;
        MAX_AMPLITUDE = pow(2, depth - 1 /* Leave a bit for sign */);
        incrementTime();
    }

    void shift(double ph) {
      phase = ph;
    }

    void start(double frequency) {
      switch(waveType) {
        case SINE:
          rps = 6.283185307179586 * frequency;
          break;
        case SQUARE:
          rps = 6.283185307179586 * frequency;
          break;
        case SAWTOOTH:
          rps = frequency;
          break;
      }
      cout << rps << endl;
    }

    double reverse(double oscillation) {
      switch(waveType) {
        case SINE:
          return asin(oscillation);
      }
    }

    double oscillate() {
        switch(waveType) {
            case SINE:
                return sine(time, rps, phase);
            case SQUARE:
                return square(time, rps, phase);
            case SAWTOOTH:
                return sawtooth(time, rps);
            case TRIANGLE:
                val = triangle(time, rps, phase, val);
                return val;
        }
        return 0.0;
    }

    void skipTo(double amplitude, bool advance) {
      double osc = oscillate();
      if(advance) incrementTime();
      stream.push_back(MAX_AMPLITUDE * amplitude * osc);
    }

    double apply(double coefficient, double offset, bool advance) {
        double amp = coefficient * time + offset;
        double osc = oscillate();
        if(advance) incrementTime();
        long res = MAX_AMPLITUDE * amp * osc;
        stream.push_back(res);
        return amp;
    }

    void pop() {
        stream.pop_back();
    }

    void incrementTime() {
        time = ++counter / sampleSize;
    }
};

class Envelope {

    public:
    double attack; //Time from key press until maximum level
    double decay; //Time from maximum level to sustain level
    double sustain; //Sustain level
    double hold; //The time to keep the sustain level
    double release; //Time from sustain until no sound
    double currentAmp;
    double holdTime;
    Envelope(double a, double d, double s, double h, double r) {
        holdTime = a + d + h;
        attack = a;
        decay = d;
        sustain = s;
        hold = h;
        release = r;
        if(r == 0.0 && sustain == 0.0)
          release = decay;
    }

    int passThrough(Oscillator &osc, double amplitude, double time = -1) {
        if(time < 0)
          time = holdTime;
        if(amplitude > 1.0 || sustain > amplitude)
            return 1;
        doAttack(osc, amplitude, time);
        if(osc.time >= time)
          return doRelease(osc, amplitude);
        doDecay(osc, amplitude, time);
        if(osc.time >= time)
          return doRelease(osc, amplitude);
        doHold(osc, amplitude, time);
        doRelease(osc, amplitude);
        return 0;
    }

    double mockAttack(Oscillator &osc, double amplitude) {
      double coef = amplitude / attack;
      return coef * osc.time;
    }

    int doAttack(Oscillator &osc, double amplitude, double time) {
        if(attack == 0.0) {
          osc.skipTo(amplitude, false);
          return 0;
        }
        double coef = amplitude / attack;
        double amp = 0;
        while(amp < amplitude && osc.time < time) {
            amp = osc.apply(coef, 0, true);
        }
        currentAmp = amp;
        osc.pop();
        return 0;
    }

    int doDecay(Oscillator &osc, double amplitude, double time) {
      if(decay == 0.0) {
        osc.skipTo(sustain, false);
        return 0;
      }
        double coef = (sustain - amplitude) / decay;
        //Passes through (attack, amplitude)
        double offset = amplitude - attack * coef;
        double amp = amplitude;
        while(amp > sustain && osc.time < time) {
            amp = osc.apply(coef, offset, true);
        }
        currentAmp = amp;
        osc.pop();
        return 0;
    }

    int doHold(Oscillator &osc, double amplitude, double time) {
        for(int i = 0; i < osc.sampleSize * hold; i++) {
            osc.apply(0, sustain, true);
            if(osc.time >= time)
              break;
        }
        if(hold == 0.0) {
          while(osc.time + release < time) {
            osc.apply(0, currentAmp, true);
          }
        }
        return 0;
    }

    int doRelease(Oscillator &osc, double amplitude) {
        if(release == 0.0) {
          osc.skipTo(0, false);
          return 0;
        }
        double coef = -(currentAmp / release);
        //Passes through (currentAmp, osc.time)
        double offset = currentAmp - osc.time * coef;
        double amp = sustain;
        while(amp > 0) {
            amp = osc.apply(coef, offset, true);
        }
        osc.pop();
        return 0;
    }
};

vector<long> mix(vector<vector<long>> oscs) {
  vector<long> res(oscs[0].size(), 0);
  for(int i = 0; i < oscs.size(); i++) {
    for(long j = 0; j < oscs[i].size(); j++) {
      res[j] += oscs[i][j];
    }
  }
  long f = oscs.size();
  for(int i = 0; i < res.size(); i++) {
    res[i] /= f;
  }
  return res;
}

vector<long> series(vector<vector<long>> oscs, bool overlap = false, vector<double> durations = {}, int sampleSize = 0, double quarter = 0.0) {
  vector<long> res;
  double duration = 0.0;
  int startIndex = 0;
  for(int i = 0; i < oscs.size(); i++) {
    if(overlap) {
      if(i == 0) {
        for(int j = 0; j < oscs[i].size(); j++) {
          res.push_back(oscs[i][j]);
        }
      } else {
        int sampleStart = sampleSize * duration;
        for(int j = 0; j < oscs[i].size(); j++) {
          int sample = sampleStart + j;
          if(sample < res.size()) {
            long avg = (res[sample] + oscs[i][j]) / 2;
            res[sample] = max(avg, res[sample]);
          } else {
             res.push_back(oscs[i][j]);
           }
        }
      }
      double quarters = 4 / durations[i];
      duration += quarters * quarter;
    } else {
      for(int j = 0; j < oscs[i].size(); j++) {
        res.push_back(oscs[i][j]);
      }
    }
  }
  /*if(overlap) {
    for(int i = 0; i < res.size(); i++) {
      res[i] = res[i] / atIndex[i];
    }
  }*/
  return res;
}

struct WAV_HEADER {
    /* RIFF Chunk Descriptor */
    uint8_t RIFF[4] = {'R', 'I', 'F', 'F'}; // RIFF Header Magic header
    uint32_t ChunkSize;                     // RIFF Chunk Size
    uint8_t WAVE[4] = {'W', 'A', 'V', 'E'}; // WAVE Header
    /* "fmt" sub-chunk */
    uint8_t fmt[4] = {'f', 'm', 't', ' '}; // FMT header
    uint32_t Subchunk1Size = 16;           // Size of the fmt chunk
    uint16_t AudioFormat = 1; // Audio format 1=PCM,6=mulaw,7=alaw,     257=IBM
                              // Mu-Law, 258=IBM A-Law, 259=ADPCM
    uint16_t NumOfChan = 1;   // Number of channels 1=Mono 2=Sterio
    uint32_t SamplesPerSec = 16000;   // Sampling Frequency in Hz
    uint32_t bytesPerSec = 16000 * 2; // bytes per second
    uint16_t blockAlign = 2;          // 2=16-bit mono, 4=16-bit stereo
    uint16_t bitDepth = 16;      // Number of bits per sample
    /* "data" sub-chunk */
    uint8_t Subchunk2ID[4] = {'d', 'a', 't', 'a'}; // "data"  string
    uint32_t Subchunk2Size;                        // Sampled data length
};

int writeWAV(vector<long> &data, int sampleSize, int bitDepth) {
    uint32_t dataSize = data.size() * (bitDepth / 8);

    WAV_HEADER header;
    header.bitDepth = bitDepth;
    header.SamplesPerSec = sampleSize;
    header.blockAlign = bitDepth / 8;
    header.bytesPerSec = header.SamplesPerSec * header.blockAlign;
    header.ChunkSize = dataSize + sizeof(WAV_HEADER) - 8;
    header.Subchunk2Size = dataSize + sizeof(WAV_HEADER) - 44;

    ofstream out("test.wav", ios::binary);
    out.write(reinterpret_cast<const char *>(&header), sizeof(header));

    for(int i = 0; i < data.size(); i++) {
      out.write(reinterpret_cast<char *>(&(data[i])), bitDepth / 8);
    }
    out.close();
    return 0;
}

vector<long> simultaneous(vector<double> notes, int sampleSize, int bitDepth, WaveType wav, Envelope env, double amplitude) {
  vector<vector<long>> res;
  for(int i = 0; i < notes.size(); i++) {
    Oscillator osc(wav, sampleSize, bitDepth);
    osc.start(notes[i]);
    env.passThrough(osc, amplitude);
    res.push_back(osc.stream);
  }
  return mix(res);
}

vector<long> successive(vector<double> notes, vector<double> durations, int sampleSize, int bitDepth, WaveType wave, Envelope env, double amplitude, double quarter) {
  double MAX_AMPLITUDE = pow(2, bitDepth - 1);
  double duration = 0.0;
  vector<long> res;
  double prevPhase = 0.0;
  double prevFreq = 0.0;
  double prevTime = 0.0;
  for(int i = 0; i < notes.size(); i++) {
    double noteDur = quarter * (4 / durations[i]);
    prevFreq = notes[i];
    if(notes[i] < 0.0) {
      long samples = sampleSize * noteDur;
      long sampleStart = sampleSize * duration;
      for(int j = 0; j < samples; j++) {
        if(j >= res.size()) res.push_back(0.0);
      }
    } else {
      long sampleStart = sampleSize * duration;
      Oscillator osc(wave, sampleSize, bitDepth);
      osc.start(notes[i]);
      /*if(wave == SINE && res.size() > sampleStart) {
        //double amp = env.mockAttack(osc, amplitude);
        for(long i = sampleStart; i < res.size(); i++) {

        }
        double phase = prevFreq * prevTime + prevPhase - osc.time * osc.rps;
        osc.shift(phase);
        prevPhase = phase;
      }*/
      env.passThrough(osc, amplitude, noteDur);
      for(int j = 0; j < osc.stream.size(); j++) {
        if(sampleStart + j >= res.size() || sampleStart + j == 0) {
          res.push_back(osc.stream[j]);
          continue;
        }
        /*long prev = res[sampleStart + j - 1];
        long nd = abs(osc.stream[j] - prev);
        long od = abs(res[sampleStart + j] - prev);
        if(nd < od) res[sampleStart + j] = osc.stream[j];*/
        //res[sampleStart + j] = min(abs(osc.stream[j] - prev), abs(res[sampleStart + j] - prev));
        res[sampleStart + j] += osc.stream[j];
        res[sampleStart + j] /= 2;
      }
    }
    duration += noteDur;
    prevTime = noteDur;
  }
  return res;
}

class Instrument {
  public:
  double amplitude;
  Envelope envelope;
  WaveType waveType;

  Instrument(double amp, Envelope &env, WaveType wav): amplitude(amp), envelope(env), waveType(wav) {
    amplitude = amp;
    envelope = env;
    waveType = wav;
  }

  vector<long> playHarmony(vector<double> notes, int sampleSize, int bitDepth) {
    return simultaneous(notes, sampleSize, bitDepth, waveType, envelope, amplitude);
  }

  vector<long> playMelody(vector<double> notes, vector<double> duration, int sampleSize, int bitDepth) {
    return successive(notes, duration, sampleSize, bitDepth, waveType, envelope, amplitude, 0.5);
  }
};

vector<long> everybreath(int sampleSize, int bitDepth) {
  Envelope env(0.01, 0.7, 0.5, 3.0, 0.0);

  vector<double> notes = {A2, E3, B3, E3, DB4, B3, E3, B3, A2, E3, B3, E3, DB4, B3, E3, B3, GB2, DB3, AB3, DB3, A3, AB3, DB3, AB3, GB2, DB3, AB3, DB3, A3, AB3, DB3, AB3, D3, A3, E4, D3, D4, A3, D3, A3, E3, B3, GB4, E3, E4, B3, E3, B3, A2, E3, B3, E3, DB4, B3, E3, B3, A2, E3, B3, E3, DB4, B3, E3, B3, A2, E3, B3, E3, DB4, B3, E3, B3, A2, E3, B3, E3, DB4, B3, E3, B3, GB2, DB3, AB3, DB3, A3, AB3, DB3, AB3, GB2, DB3, AB3, DB3, A3, AB3, DB3, AB3, D3, A3, E4, D3, D4, A3, D3, A3, E3, B3, GB4, E3, E4, B3, E3, B3};
  vector<double> durations = {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
  return successive(notes, durations, sampleSize, bitDepth, SQUARE, env, 0.5, 0.51282);
}

vector<long> crazytrain(int sampleSize, int bitDepth) {
  Envelope env(0.01, 0.7, 0.0, 0.0, 0.0);

  vector<double> notes = {GB2, GB2, DB3, GB2, D3, GB2, DB3, GB2, B2, A2, AB2, A2, B2, A2, AB2, E2, GB2, GB2, DB3, GB2, D3, GB2, DB3, GB2, B2, A2, AB2, A2, B2, A2, AB2, E2, GB2, GB2, DB3, GB2, D3, GB2, DB3, GB2, B2, A2, AB2, A2, B2, A2, AB2, E2};
  vector<double> durations = {8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8};
  return successive(notes, durations, sampleSize, bitDepth, SQUARE, env, 0.5, 0.434782608695);
}

vector<long> cantstop(int sampleSize, int bitDepth) {
  Envelope env(0.01, 1.2, 0.0, 0.0, 0.0);

  vector<double> notes = {D3, E3, D3, E3, D3, E3, D3, E3, D3, E3, 0.0, D3, E3, D3, E3, D3, E3, D3, E3, D3, E3, 0.0, D3, E3, D3, E3, D3, E3, D3, E3, D3, E3, 0.0, D3, E3, D3, E3, D3, E3, D3, E3, D3, E3, 0.0, D3, E3, D3, E3, D3, E3, D3, E3, D3, E3, 0.0, D3, E3, D3, E3, D3, E3, D3, E3, D3, E3, 0.0, D3, E3, D3, E3, D3, E3, D3, E3, D3, E3, 0.0, D3, E3, D3, E3, D3, E3, D3, E3, D3, E3, 0.0};
  vector<double> durations = {16, 8, 16, 8, 16, 8, 16, 8, 16, 16, 8, 16, 8, 16, 8, 16, 8, 16, 8, 16, 16, 8, 16, 8, 16, 8, 16, 8, 16, 8, 16, 16, 8, 16, 8, 16, 8, 16, 8, 16, 8, 16, 16, 8, 16, 8, 16, 8, 16, 8, 16, 8, 16, 16, 8, 16, 8, 16, 8, 16, 8, 16, 8, 16, 16, 8, 16, 8, 16, 8, 16, 8, 16, 8, 16, 16, 8, 16, 8, 16, 8, 16, 8, 16, 8, 16, 16, 8};
  return successive(notes, durations, sampleSize, bitDepth, SQUARE, env, 0.5, 0.674);
}

vector<long> wearenumberone(int sampleSize, int bitDepth) {
  Envelope env(0.1, 0.3, 0.9, 0.0, 0.5);

  vector<double> notes = {F3, 0.0, C4, B3, C4, B3, C4, B3, C4, AB3, F3, 0.0, F3, 0.0, AB3, 0.0, C4, 0.0, DB4, 0.0, AB3, 0.0, DB4, 0.0, EB4, 0.0, C4, DB4, C4, DB4, C4, 0.0, 0.0, F3, 0.0, C4, B3, C4, B3, C4, B3, C4, AB3, F3, 0.0, F3, 0.0, AB3, 0.0, C4, 0.0, DB4, 0.0, AB3, 0.0, DB4, 0.0, EB4, 0.0, C4, DB4, C4, DB4, C4, 0.0, 0.0, F3, 0.0, C4, B3, C4, B3, C4, B3, C4, AB3, F3, 0.0, F3, 0.0, AB3, 0.0, C4, 0.0, DB4, 0.0, AB3, 0.0, DB4, 0.0, EB4, 0.0, C4, DB4, C4, DB4, C4, 0.0, 0.0};
  vector<double> durations = {4, 8, 8, 16, 16, 16, 16, 8, 8, 4, 4, 8, 16, 16, 16, 16, 16, 16, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 4, 8, 8, 16, 16, 16, 16, 8, 8, 4, 4, 8, 16, 16, 16, 16, 16, 16, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8, 4, 8, 8, 16, 16, 16, 16, 8, 8, 4, 4, 8, 16, 16, 16, 16, 16, 16, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 4, 8};
  return successive(notes, durations, sampleSize, bitDepth, SQUARE, env, 1.0, 0.4);
}

vector<long> ageofwar(int sampleSize, int bitDepth) {
  //150 BPM
  Envelope env(0.07, 1.8, 0.3, 0.0, 0.01);

  /*Oscillator osc(SQUARE, sampleSize, bitDepth);
  osc.start(D5);
  quarter.passThrough(osc, 0.4);
  res = osc.stream;
  vector<double> fc = {B4, GB4};
  Envelope third(0.0, 0.0, 0.4, 0.133333, 0.0);
  res = series({res, simultaneous(fc, sampleSize, bitDepth, SINE, third, 0.4)});*/
  //vector<double> notes = {B2, D3, E3, GB3, A3, B3, DB4, E4, GB4, A4, B4, D5, B4, A4, B4, GB4, A4, B4, E4, A4, E4, DB3, D3, GB4, D3, E3, D3, DB3, D3, DB3, D3, A4, B3, A3, E4, D4, B4, A4, B4, GB4, A4, B4, E4, A4, E4, DB3, D3, GB4, D3, E3, D3, DB3, D3, DB3, D3, A4, B3, A3, E3, B4, D3, E3, GB3, A3, B3, DB4, E4, GB4, A4, B4, D5};
  vector<double> notes = {B3, D4, E4, GB4, A4, B4, DB5, E5, GB5, A5, B5, D6, B3, A3, B3, GB3, A3, B3, E3, A3, E3, DB4, D4, GB3, D4, E4, D4, DB4, D4, DB4, D4, A3, B4, A4, E3, D3, B3, A3, B3, GB3, A3, B3, E3, A3, E3, DB4, D4, GB3, D4, E4, D4, DB4, D4, DB4, D4, A3, B4, A4, E4, B3, D4, E4, GB4, A4, B4, DB5, E5, GB5, A5, B5, D6};
  vector<double> durations = {12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 4, 4, 6, 12, 4, 4, 4, 6, 12, 4, 4, 4, 4, 6, 12, 4, 4, 6, 12, 4, 1, 1, 1, 1, 4, 4, 6, 12, 4, 4, 4, 6, 12, 4, 4, 4, 4, 6, 12, 4, 4, 6, 12, 4, 1, 1, 1, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12};
  /*vector<double> notes = {D4, B4, B4, -1, D4, D3, D3, -1, E4};
  vector<double> durations = {4, 6, 12, 2, 4, 6, 12, 4, 4};*/
  vector<long> res = successive(notes, durations, sampleSize, bitDepth, SQUARE, env, 1, 0.4);
  return res;
}

vector<long> sevennationarmy(int sampleSize, int bitDepth) {
  //124 BPM
}

vector<long> nokia(int sampleSize, int bitDepth) {
    //120 BPM
    Envelope env(0.05, 0.3, 0.0, 0.0, 0.1);
    Instrument ins(0.5, env, SQUARE);

    Oscillator osc(SINE, sampleSize, bitDepth);

    vector<double> notes = {E4, D4, GB3, AB3, DB4, B3, D3, E3, B3, A3, DB3, E3, A3};
    vector<double> durations = {16, 16, 8, 8, 16, 16, 8 , 8, 16, 16, 8, 8, 1};
    vector<long> res = successive(notes, durations, sampleSize, bitDepth, SQUARE, env, 0.5, 0.5);
    return res;
}

int main() {
    int sampleSize = 44100;
    int bitDepth = 24;

    vector<long> res = everybreath(sampleSize, bitDepth);

    writeWAV(res, sampleSize, bitDepth);
    return 0;
}
