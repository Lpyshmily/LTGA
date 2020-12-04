clear;clc;
raw_data = dlmread('trajectory.txt');
data = raw_data(:, 2:7);

r = data(:, 1:3);

figure(1)
plot3(r(:,1), r(:,2), r(:,3));
% axis equal
figure(2)
plot(r(:,1), r(:,2))
