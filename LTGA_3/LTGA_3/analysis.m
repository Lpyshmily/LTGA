clear;clc;
data = dlmread('info_1606383602.txt');
% [min_value, min_index] = min(data(:,4));
% xbest = data(min_index, 1:3);
% fbest = data(min_index, 4);
valid_index = data(:,4) < 0;
valid_data = data(valid_index, :);

figure(1)
scatter3(valid_data(:,1), valid_data(:,2), valid_data(:,3), 30, valid_data(:,4), 'filled');
colorbar

figure(2)

data_1_index = valid_data(:,4) < -1.56e4;
data_1 = valid_data(data_1_index, :);
scatter3(data_1(:,1), data_1(:,2), data_1(:,3), 30, 'filled', 'r');
hold on

xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

data_2_index = (valid_data(:,4) < -1.54e4) & (valid_data(:,4) > -1.56e4);
data_2 = valid_data(data_2_index, :);
scatter3(data_2(:,1), data_2(:,2), data_2(:,3), 30, 'filled', 'b');
hold on

data_3_index = (valid_data(:,4) < -1.52e4) & (valid_data(:,4) > -1.54e4);
data_3 = valid_data(data_3_index, :);
scatter3(data_3(:,1), data_3(:,2), data_3(:,3), 30, 'filled', 'g');
hold on

data_4_index = (valid_data(:,4) < -1.50e4) & (valid_data(:,4) > -1.52e4);
data_4 = valid_data(data_4_index, :);
scatter3(data_4(:,1), data_4(:,2), data_4(:,3), 30, 'filled');
hold on

data_5_index = valid_data(:,4) > -1.50e4;
data_5 = valid_data(data_5_index, :);
scatter3(data_5(:,1), data_5(:,2), data_5(:,3), 30, 'filled');