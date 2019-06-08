Fs = 4e3; Ts=1/Fs; tones = [];
t = 0:Ts:0.5-Ts;
ver = [697 770 852 941];
hor = [1209 1336 1477];
for k = 1:length(ver)
    for l = 1:length(hor)
        tone = sum(sin(2*pi*[ver(k);hor(l)].*t))';
        tones = [tones;tone;zeros(size(tone))];
    end
end
soundsc(tones,Fs); % To hear
S = timetable(seconds(0:length(tones)-1)'/Fs,tones);
total_time = seconds(0:numel(tones)-1)'/Fs;
pspectrum(tones,total_time,'spectrogram','TimeResolution',0.2,'OverlapPercent',99,'Leakage',0.85)