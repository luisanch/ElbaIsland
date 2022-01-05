clear variables
close all
clc

load('received_signals.mat')

%% inputs
h = 150; %distance to ocean floor fromthe surface
hr = 15 : 15 : h - 15; %z positions of the transceivers
receivers = [(zeros(size(hr)))', hr'];  %array with x,z positions of the transceivers
n_s = 8; %amount of computed sources
sampleTime = 8.3361e-04; 
c = 1500; %speed of sound
resolution = 10; %resolution in meters for back prpagation in selected area
yRange = (0 : resolution : h);  %vertical range of selected area
xRange = (0 : resolution : 1500);  %horizontal range of selected area

%% inputs for simulated signal 
source = [50 50]; %source loccation
a = 10; %source signal amplitude
freq = 10; %in hz 
totalTime = 1.6; % in s

%% Simulate Data
%s_t = emit(totalTime, sampleTime, @mySinc, freq, a);
%n_gg_t= fitArrays(getMultipleRecordings(source, h, n_s, receivers, sampleTime, s_t, c)); 
%sr_t = backPropagate(n_gg_t, receivers, source, sampleTime, n_s, c, h);

%% Plots for simulated data
% plot(s_t(:,1), s_t(:,2))
% title("original")
% figure()
% plot(sr_t(:,1), sr_t(:,2))
% title("back Propagated") 

%% Get Real Data 
n_gg_t= greenToRecorded(green, sampleTime);

%% Plots 
figure()
plotAll(n_gg_t, false)
title("recorded")
grid = getMultipleBackPropagation(n_gg_t, receivers, sampleTime, n_s, c, h, yRange, xRange);
figure()
image(grid,'CDataMapping','scaled') 
xticklabels((xticks*resolution)) 
yticklabels((yticks*resolution))
colorbar

%% back propagate in multiple points
function grid = getMultipleBackPropagation(n_gg_t, receivers, sampleTime, n_s, c, h, yRange, xRange)
grid = zeros([size(yRange,2) size(xRange,2)]);
for x = 1 : 1 : size(xRange,2)
    for y = 1 : 1 : size(yRange,2)
        sr_t = backPropagate(n_gg_t, receivers, [xRange(x) yRange(y)], sampleTime, n_s, c, h);
        grid(y,x) = max(sr_t(:,2));
    end
end
end

%% back propagation
function sr_t = backPropagate(n_gg_t, receivers, atLocation, sampleTime, n_s, c, h)
n_ss_t = {};
for i =  1 : 1 : size(n_gg_t, 2)
    reverse_source = receivers(i,:);
    reverse_receiver = atLocation;
    reverse_s_t = [n_gg_t{i}(:,1), flipud(n_gg_t{i}(:,2))];
    n_ss_t{i} = sumArrays(fitArrays(recordMultSources(reverse_source, h, n_s, reverse_receiver, sampleTime, reverse_s_t, c)));
end
sr_t = normalize(sumArrays(fitArrays(n_ss_t))); %not sure about this
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
    gg_t{j} = record(s_t, source, receiver, c, sampleTime, j-1, h);
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

g_t(:,1) =  g_t(:,1) - g_t(1,1); %adjust delay
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
e = 1;
%%commented this out for the time being since e = +1
% if ~mod(j,2) % j is even
%     e = (-1)^(j/2);
% else         % j is odd
%     e = (-1)^((j-1)/2);
% end
end

%% function generation utils
function n_gg_t = greenToRecorded(green, sampleTime)
n_gg_t = {};
for i = 1 : 1 : size(green,1)
    n_gg_t{i} = [(0:1:size(green,2)-1)',green(i,:)'];
end
end

function sinc = mySinc(t)
if t == 0
    sinc = 1;
else
    sinc = sin(t)/(t);
end
end

function s_t = emit(totalTime, sampleTime, func, freq, a)
freqRad = pi * freq; %half
sampleSpace = 0 : sampleTime : totalTime/2;
sampleSpace = [flip(-sampleSpace) sampleSpace];
s_t = zeros(size(sampleSpace,2), 2);

for i = 1 : 1 : size(sampleSpace,2)
    s_t(i, 1) =  sampleSpace(i) + sampleSpace(size(sampleSpace,2));
    s_t(i, 2) = propagateWave(func,freqRad, sampleSpace(i), a);
end
end

% function s_t = emit(totalTime, sampleTime, func, freq, a)
% freqRad = 2 * pi * freq;
% sampleSpace = 0 : sampleTime : totalTime;
% s_t = zeros(totalTime / sampleTime, 2);
%
% for i = 1 : 1 : totalTime / sampleTime
%     s_t(i, 1) = sampleSpace(i);
%     s_t(i, 2) = propagateWave(func,freqRad, sampleSpace(i), a);
% end
% end

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
function y = plotAll(cellArray, singleWindow)
legendLabel = {}; 
for i = 1 : 1 : size(cellArray,2) 
    if singleWindow
    legendLabel{i} = num2str(i);
    plot(cellArray{i}(:,1), cellArray{i}(:,2));
    hold on
    else
    subplot(size(cellArray,2), 1, i)
    plot(cellArray{i}(:,1), cellArray{i}(:,2));
    legend({num2str(i)},'Location','southwest') 
    end
end
if singleWindow
legend(legendLabel,'Location','southwest')
end
end