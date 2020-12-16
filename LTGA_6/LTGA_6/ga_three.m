clear;clc;
% data = dlmread('ga_three_obj_nlopt_list.txt'); % [11*11*21,4]
% sortedData = zeros(size(data));
% for i=1:11*11
%     sortedData(21*i-20:21*i, :) = sortrows(data(21*i-20:21*i, :),3);
% end
% 
% valid_index = sortedData(:,4) < 0;
% valid_data = sortedData(valid_index, :);

data = dlmread('ga_three_obj_nlopt_denselist.txt'); % [11*11*11,4]
sortedData = zeros(size(data));
for i=1:11*11
    sortedData(11*i-10:11*i, :) = sortrows(data(11*i-10:11*i, :),3);
end

valid_index = sortedData(:,4) < 0;
valid_data = sortedData(valid_index, :);


figure(1)
scatter3(valid_data(:,1), valid_data(:,2), valid_data(:,3), 30, valid_data(:,4), 'filled');
colorbar

figure(2)

data_1_index = valid_data(:,4) < -1.56e4;
data_1 = valid_data(data_1_index, :);
scatter3(data_1(:,1), data_1(:,2), data_1(:,3), 30, 'filled', 'r');
hold on

xlim([0.3,0.35]);
ylim([0,1]);
zlim([0,0.5]);

data_2_index = (valid_data(:,4) < -1.54e4) & (valid_data(:,4) > -1.56e4);
data_2 = valid_data(data_2_index, :);
scatter3(data_2(:,1), data_2(:,2), data_2(:,3), 30, 'filled', 'g');
hold on

data_3_index = (valid_data(:,4) < -1.52e4) & (valid_data(:,4) > -1.54e4);
data_3 = valid_data(data_3_index, :);
scatter3(data_3(:,1), data_3(:,2), data_3(:,3), 30, 'filled', 'b');
hold on

data_4_index = (valid_data(:,4) < -1.50e4) & (valid_data(:,4) > -1.52e4);
data_4 = valid_data(data_4_index, :);
scatter3(data_4(:,1), data_4(:,2), data_4(:,3), 30, 'filled', 'm');
hold on

data_5_index = valid_data(:,4) > -1.50e4;
data_5 = valid_data(data_5_index, :);
scatter3(data_5(:,1), data_5(:,2), data_5(:,3), 30, 'filled', 'k');

figure(3)
valid_data(:,4) = -valid_data(:,4);
subplot(3,2,1);
data30 = valid_data(valid_data(:,1)==0.30, :);
scatter(data30(:,2), data30(:,3), 30, data30(:,4), 'filled');

subplot(3,2,2);
data31 = valid_data(valid_data(:,1)==0.31, :);
scatter(data31(:,2), data31(:,3), 30, data31(:,4), 'filled');

subplot(3,2,3);
data32 = valid_data(valid_data(:,1)==0.32, :);
scatter(data32(:,2), data32(:,3), 30, data32(:,4), 'filled');

subplot(3,2,4);
data33 = valid_data(valid_data(:,1)==0.33, :);
scatter(data33(:,2), data33(:,3), 30, data33(:,4), 'filled');

subplot(3,2,5);
data34 = valid_data(valid_data(:,1)==0.34, :);
scatter(data34(:,2), data34(:,3), 30, data34(:,4), 'filled');

subplot(3,2,6);
data35 = valid_data(valid_data(:,1)==0.35, :);
scatter(data35(:,2), data35(:,3), 30, data35(:,4), 'filled');

PSO_data = dlmread('ga_three_PSO.txt');
figure(4)
plot(PSO_data);
