% MPPI for Lane-Following Control of Vehicle
clc;
clear;

% Parameters
L = 2.5; % wheelbase (m)
Lr = 1.25; % distance from rear axle to CoG (m)
Ts = 0.01; % sampling time (s)
horizon = 10; % control horizon
num_samples = 10; % number of samples for MPPI
lambda = 1.0; % temperature parameter for MPPI
sigma = [2; 2]; % standard deviation for sampling control inputs

% Define the reference trajectory (straight line for simplicity)
refTrajectory.x = 0:0.1:100; % x-coordinates
%refTrajectory.y = 2 * ones(size(refTrajectory.x)); % y-coordinates (lane center at y=2m)
refTrajectory.y = sqrt(refTrajectory.x)+2; % y-coordinates (lane center at y=2m)
refTrajectory.v = 10*ones(length(refTrajectory.x),1); % target_velocity(m/s)

goal=[10,10];

% Initial conditions
x0 = 0; % initial x position (m)
y0 = 2; % initial y position (m)
theta0 = 0; % initial orientation (rad)
v0 = 10; % initial velocity (m/s)

% Initial state and control input
state = [x0; y0; theta0;v0];
states = zeros(4, length(refTrajectory.x));
states(:, 1) = state;
controls = [0; v0]; % initial control input: [delta; velocity]

% Store predicted paths
predicted_paths = zeros(4, horizon, length(refTrajectory.x) - 1);

% Main simulation loop
for k = 1:length(refTrajectory.x) - 1
   

    % 本来はautoware-practiceのPurepursuitみたいな感じで最近傍の距離から目標経路マップを紹介して目標位置と目標車速を出さないといけない
    
    % Generate control samples
    control_samples = repmat(controls, 1, num_samples) + sigma .* randn(2, num_samples);
    
    % Keep control signals within bounds
    control_samples(1, :) = max(min(control_samples(1, :), pi/4), -pi/4); % steering angle between -45 to 45 degrees
    control_samples(2, :) = max(min(control_samples(2, :), 20), 0); % velocity between 0 to 20 m/s

    costs = zeros(1, num_samples);
    
    % Evaluate each sample
    for i = 1:num_samples
        sample_state = state;
        cost = 0;
        valid_sample = true;
        for t = 1:horizon
            % Apply control input to the vehicle model
            u = control_samples(:, i);
            sample_state = kinematicBicycleModel(sample_state, u, Ts, L, Lr);
            
            % Check for NaN or Inf in the sample state
            if any(isnan(sample_state)) || any(isinf(sample_state))
                valid_sample = false;
                break;
            end

            min_distance=10000000; % way point最近傍点の初期化
            for j = 1:length(refTrajectory.x) - 1
                distance=sqrt((state(1)-refTrajectory.x(j))^2+(state(2)-refTrajectory.y(j))^2); % ユークリッド距離
                if min_distance>distance
                    closest_way_point=j;
                    min_distance=distance;
                end
            end
            % Current reference position
            ref = [refTrajectory.x(closest_way_point); refTrajectory.y(closest_way_point); refTrajectory.v(closest_way_point)];
            
            % Compute cost (e.g., distance to reference)
            cost = cost + norm(sample_state(1:2) - ref(1:2), 2)+1*norm(sample_state(4)-refTrajectory.v(closest_way_point),1)+ 0.1*norm(u, 2); 
            %refの取り方が不適なのでバグる。目標経路と自車位置の最近傍点を計算し、曲率と道なり距離から将来目標位置を出すべき
%            if (sample_state(1)>60)&&(sample_state(1)<80)&&(sample_state(2)<10)&&(sample_state(2)>-10)
%                cost=cost+1000;
%            end


        end
        cost=cost+100*norm(sample_state(1:2) - ref(1:2), 2)+10*norm(sample_state(4)-refTrajectory.v(closest_way_point),1);

        if valid_sample
            costs(i) = cost;
        else
            costs(i) = inf;
        end
    end
    
    % Check for NaN or Inf in costs
    if any(isnan(costs)) || all(isinf(costs))
        error('NaN or Inf detected in costs. Check model dynamics and cost calculation.');
    end
    
    % Update control input based on costs
    weights = exp(-costs / lambda);
    weights = weights / sum(weights);
    controls = sum(control_samples .* weights, 2);
    
    % Predict the path using the best control input
    predicted_state = state;
    for t = 1:horizon
        predicted_state = kinematicBicycleModel(predicted_state, controls, Ts, L, Lr);
        predicted_paths(:, t, k) = predicted_state;
    end
    
    % Apply the control input to the vehicle model
    state = kinematicBicycleModel(state, controls, Ts, L, Lr);
    states(:, k + 1) = state;
end

% Plot results
figure;
plot(refTrajectory.x, refTrajectory.y, 'r--', 'LineWidth', 2);
hold on;
plot(states(1, :), states(2, :), 'b', 'LineWidth', 2);

% Plot predicted paths
for k = 1:length(refTrajectory.x) - 1
    plot(predicted_paths(1, :, k), predicted_paths(2, :, k), 'g--');
end

xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Lane-Following Control using MPPI with Predicted Paths');
legend('Reference Path', 'Actual Path', 'Predicted Path');
grid on;

% Kinematic Bicycle Model Function
function state = kinematicBicycleModel(state, u, Ts, L, Lr)
    % Unpack state and control input
    x = state(1);
    y = state(2);
    theta = state(3);
    v = u(2);
    delta = u(1);
    
    % Kinematic bicycle model equations
    beta = atan(Lr * tan(delta) / L);
    if isnan(beta) || isinf(beta)
        beta = 0;
    end
    x = x + Ts * v * cos(theta + beta);
    y = y + Ts * v * sin(theta + beta);
    theta = theta + Ts * v * sin(beta) / L;
    
    % Update state
    state = [x; y; theta;v];
end
