data = dlmread('info_triple.txt');
valid_index = data(:,4) < 0;
valid_data = data(valid_index, :);
scatter3(valid_data(:,1), valid_data(:,2), valid_data(:,3), 50, valid_data(:,4));
colorbar