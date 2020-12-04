clear;clc;
data = dlmread('amplitude.txt');
data(1001:2000, 1) = data(1001:2000, 1) + 14.682744;
plot(data(:,1), data(:,2),'.-');