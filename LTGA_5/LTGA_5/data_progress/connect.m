clear; clc;
trajectory_data = dlmread('trajectory.txt');
trajectory_data(888:end, 1) = trajectory_data(888:end, 1) + 14.6827436227285;
amplitude_data = dlmread('amplitude.txt');
amplitude_data(1001:2000, 1) = amplitude_data(1001:2000, 1) + 14.682744;

timeList = trajectory_data(:,1);
time_length = size(trajectory_data, 1);
r = trajectory_data(:, 2:4);

% 寻找推力大小跳变的位置
jump_index = 1;
current_state = amplitude_data(jump_index,2);
while(1)
    jump_index = jump_index + 1;
    if (amplitude_data(jump_index,2)) ~= current_state
        break;
    end
end
jump_time = amplitude_data(jump_index,1);

start_index = 1;
stop_index = 0;
flag = 1; % 先满推设置成1，先滑行设置成0
figure(1)
for i=1:time_length
    
    if timeList(i)>jump_time
        % stop_index = i-1;
        stop_index = i;
        if flag==0
            plot(r(start_index:stop_index,1), r(start_index:stop_index,2), 'b');
            hold on
            flag = 1;
        else
            plot(r(start_index:stop_index,1), r(start_index:stop_index,2), 'r');
            hold on
            flag = 0;
        end
        start_index = i;
        % 寻找下一个跳变点
        current_state = amplitude_data(jump_index,2);
        while(1)
            jump_index = jump_index + 1;
            if (jump_index<=size(amplitude_data,1) && amplitude_data(jump_index,2)) ~= current_state
                break;
            end
        end
        jump_time = amplitude_data(jump_index-1,1);
    end
end
% 最后一段轨迹
stop_index = size(trajectory_data,1);
if flag==0
    plot(r(start_index:stop_index,1), r(start_index:stop_index,2), 'b');
else
    plot(r(start_index:stop_index,1), r(start_index:stop_index,2), 'r');
end