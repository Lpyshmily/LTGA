clear;clc;
data = dlmread('ga_single1_denselist.txt');
plot(data(:,1), -data(:,2));