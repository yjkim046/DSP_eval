% Rotor representation with Kalman Filter

Tsc = 0.1; wb = 0.5;
q = 1; Q = q*eye(3);
r = 1; R = r*eye(2);

% Initialization
x = [1 0 0.1]'; % [cos(th) sin(th) w]
u = [0 0]';     % No external force
P_h = eye(3);   % Covariance matrix

x_k = zeros(3,1000);
QF = 10;
for k = 1:size(x_k,2)
    % Projection
    x_p = [x(1)*cos(wb*Tsc*x(3))-x(2)*sin(wb*Tsc*x(3));
        x(1)*sin(wb*Tsc*x(3))+x(2)*cos(wb*Tsc*x(3));
        x(3)];
    x_p = round(x_p*2^QF)/2^QF;
    F = [cos(wb*Tsc*x(3)) -sin(wb*Tsc*x(3)) -x(2)*wb*Tsc*cos(wb*Tsc*x(3));
        sin(wb*Tsc*x(3)) cos(wb*Tsc*x(3)) -x(2)*wb*Tsc*sin(wb*Tsc*x(3));
        0 0 1];
    F = round(F*2^QF)/2^QF;
    P_p = F*P_h*F' + Q;
    P_p = round(P_p*2^QF)/2^QF;

    % Kalman Gain Computation
    H = [1 0 0;0 1 0];
    K = P_p*H'/(H*P_p*H' + R);

    % Estimation and update
    y_p = H*x_p;
    x = x_p; % + K*(z - y_p); There is no noisy measurement
    P_h = P_p - K*H*P_p;
    P_h = round(P_h*2^QF)/2^QF;
    x_k(:,k) = x;
end
close all;
figure; plot(1:size(x_k,2),x_k(1,:),1:size(x_k,2),x_k(2,:)); grid on;
xlabel('k'); ylabel('value'); legend('x(1)','x(2)');