function kalman_plot(x, y, A, Q, C, R)
    sx = zeros(1, length(x));
    sy = zeros(1, length(y));
    sx(1) = x(1);
    sy(1) = y(1);
    state = zeros(size(A,1),1);
    state(1) = x(1); 
    state(2) = y(1);
    covariance = eye(size(A,1));

    for i=2:length(x)
        [state, covariance] = kalman_update(A,C,Q,R,[x(i),y(i)]', state, covariance);
        sx(i) = state(1);
        sy(i) = state(2);
    end
    
    plot(x, y, 'o-r'); 
    hold on;
    plot(sx, sy, 'o-b')