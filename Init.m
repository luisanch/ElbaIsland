clear variables
close all
clc

source = [0 -100];
receivers = [100, -50; 100, -100; 100, -150 ; 100, -200];
h = 200;
n_s = 50; %amount of computed sources
totalTime = 1.5; % in s
sampleTime = 0.01;
a = 10;
freq = 1; %in hz
c = 1500;
s_t = emit(totalTime, sampleTime, @sin, freq, a);
n_gg_t= fitArrays(getMultipleRecordings(source, h, n_s, receivers, sampleTime, s_t, c));
sr_t = getBackPropagation(n_gg_t, receivers, source, sampleTime, n_s, c, h);
plot(sr_t(:,1), sr_t(:,2))
figure()
plot(s_t(:,1), s_t(:,2))
figure()
plotAll(n_gg_t)

%% back propagation
function sr_t = getBackPropagation(n_gg_t, receivers, source, sampleTime, n_s, c, h)
n_ss_t = {};
disp(size(n_gg_t, 2))
for i =  1 : 1 : size(n_gg_t, 2)
    reverse_source = receivers(i,:);
    reverse_h = h - receivers(i,2);
    reverse_receiver = source;
    reverse_s_t = n_gg_t{i};
    n_ss_t{i} = sumArrays(fitArrays(recordMultSources(reverse_source, reverse_h, n_s, reverse_receiver, sampleTime, reverse_s_t, c)));
end
sr_t =  sumArrays(fitArrays(n_ss_t));
end

%% record signal utils
function n_gg_t = getMultipleRecordings(source, h, n_s, receivers, sampleTime, s_t, c)
for i =  1 : 1 : size(receivers, 1)
    n_gg_t{i} = recordMultSources(source, h, n_s, receivers(i,:), sampleTime, s_t, c);
    fitted = fitArrays(n_gg_t{i});
    n_gg_t{i} = sumArrays(fitted);
end
end

function gg_t = recordMultSources(source, h, n_s, receiver, sampleTime, s_t, c)
for j = 1 : 1 : n_s
    gg_t{j} = record(s_t, source, receiver, c, sampleTime, j, h);
end
end

function g_t = record(s_t, source, receiver, c, sampleTime, j, h)
x = receiver(1);
x_s = source(1);
z = receiver(2);
z_s = source(2);

r = getR(x, x_s, z, z_s, j, h);
e = getEpsilon(j);

delay = r / c;
delayedSampleSpace = -(sampleTime : sampleTime : delay); %crete delay in time axis
delayValue = zeros(1, size(delayedSampleSpace,2)); %populate with zeros before signal arrives

g_t_0 = [flip(delayedSampleSpace); delayValue];
s_t_0 = [s_t(:, 1), s_t(:, 2)*(e/ r)];
g_t = [g_t_0'; s_t_0];

g_t(:,1) =  g_t(:,1) -  g_t(1,1); %adjust delay
end

%% green params computation utils
function r = getR(x,x_s,z,z_s, j, h)
if ~mod(j,2) % j is even
    r = ((x-x_s)^2 + (z - z_s + (j*h))^2)^(1/2);
else         % j is odd
    r = ((x-x_s)^2 + (z+z_s - ((j + 1)*h))^2)^(1/2);
end
end

function e = getEpsilon(j)
if ~mod(j,2) % j is even
    e = -1^(j/2);
else         % j is odd
    e = -1^((j-1)/2);
end
end

%% function generation utils
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

%% array manipulation utils
function sum = sumArrays(cellArray)
sum = cellArray{1};
for i = 2 : 1 : size(cellArray,2)
    sum(:,2) = sum(:,2) + cellArray{i}(:,2);
end
end

function fittedArr = fitArrays(cellArray)
lengths = zeros(size(cellArray,2),1);

%find longest array
for i = 1 : 1 : size(cellArray,2)
    lengths(i) = size(cellArray{i},1);
end
[maxLength, maxlengthIndex] = max(lengths);

for i = 1 : 1 : size(cellArray,2)
    fittedArr{i} = zeros(maxLength,2);
    fittedArr{i}(:,1) = cellArray{maxlengthIndex}(:,1);
    fittedArr{i}(1:size(cellArray{i},1),2) = cellArray{i}(:,2);
end
end

%% plot Utils
function y = plotAll(cellArray)
for i = 1 : 1 : size(cellArray,2)
    plot(cellArray{i}(:,1), cellArray{i}(:,2));
    hold on
end
end