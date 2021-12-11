clear variables
close all

sources = [0 -100];
receivers = [100 -100];
s_t = []; %table for the source
totalTime = 1.5; % in s
sampleTime = 0.01;
a = 10;
freq = 1; %in hz
c = 1500;

s_t = emit(totalTime, sampleTime, @sin, freq, a);
plot(s_t(:,1), s_t(:,2));
figure()
g_t = record(s_t, sources, receivers, c, sampleTime);
plot(g_t(:,1), g_t(:,2));
function g_t = record(s_t, source, receiver, c, sampleTime)
x = receiver(1);
x_s = source(1);
z = receiver(2);
z_s = source(2);

r = getR(x,x_s,z,z_s);
delay = r / c;
delayedSampleSpace = -(sampleTime : sampleTime : delay); %crete delay in time axis
delayValue = zeros(1, size(delayedSampleSpace,2)); %populate with zeros before signal arrives
g_t_0 = [flip(delayedSampleSpace); delayValue]; 
g_t = [g_t_0'; (s_t) / r]; 
g_t(:,1) =  g_t(:,1) -  g_t(1,1); %adjust delay
end

function r = getR(x,x_s,z,z_s)
r = ((x-x_s)^2 + (z-z_s)^2)^(1/2);
end

function s_t = emit(totalTime, sampleTime, func, freq, a)
freqRad = 2 * pi * freq;
sampleSpace = 0 : sampleTime : totalTime;
s_t = zeros(totalTime / sampleTime, 2);

for i = 1 : 1 : totalTime / sampleTime
    s_t(i, 1) = sampleSpace(i);
    s_t(i, 2) = propagateWave(func,freqRad, sampleSpace(i), a);
end

end

function y = propagateWave(func,freq, t, a)
y = a * func(freq * t);
end