function [state, location] = ms_initialize(I, region, varargin)
    %colour spaces
    %I = rgb2ycbcr(I);
    %I = rgb2hsv(I);
    %I = rgb2lab(I);
    
    %parameters
    bins = 16;
    sigma = 1;
    [height, width, ~] = size(I);

    % If the provided region is a polygon ...
    if numel(region) > 4
        x1 = round(min(region(1:2:end)));
        x2 = round(max(region(1:2:end)));
        y1 = round(min(region(2:2:end)));
        y2 = round(max(region(2:2:end)));
        region = round([x1, y1, x2 - x1, y2 - y1]);
    else
        region = round([round(region(1)), round(region(2)), ... 
            round(region(1) + region(3)) - round(region(1)), ...
            round(region(2) + region(4)) - round(region(2))]);
    end

    %extract template from image
    x1 = max(1, region(1));
    y1 = max(1, region(2));
    x2 = min(width-2, region(1) + region(3) - 1);
    y2 = min(height-2, region(2) + region(4) - 1);

    %template = I((y1:y2)+1, (x1:x2)+1, :);
    w = x2 - x1 + 1;
    h = y2 - y1 + 1;
    
    %get epanechnik kernel
    eK = create_epanechnik_kernel(w, h, sigma);
    w = size(eK,2);
    h = size(eK,1);
   
    %extract template using eK dimensions
    x = (x1):(x1 + w-1);
    y = (y1):(y1 + h-1);
    template = I(int32(y), int32(x), :);
    
    %get weighted histogram from template
    Q = extract_histogram(template, bins, eK);
    Q = Q/sum(Q(:));
    
    %background hist extraction
    bckgborder = 10;
    bckg = get_patch(I, [x1 + x2 + 1, y1 + y2 + 1] / 2, 1, [w, h]+2*bckgborder);
    [bh, bw, ~] = size(bckg);
    bckgWeights = ones(bh, bw);
    bckgWeights = outOfBorderToZero(bckgWeights, [x1 + x2 + 1, y1 + y2 + 1] / 2, bw, bh, width, height);
    bckgWeights(bckgborder:bh-bckgborder, bckgborder:bw-bckgborder, :) = 0;
    B = extract_histogram(bckg, bins, bckgWeights);
    
    %weigh so the smallest non zero value is the largest weight
    m=min(B(B>0));
    B = min((m ./ B), 1);
    B = B/sum(B(:));
    
    %coordinates for ms
    [X, Y] = meshgrid(((-w+1)/2):((w-1)/2), ((-h+1)/2):((h-1)/2)); 
    
    %create struct
    state = struct('template', template, 'size', [w, h]);
    state.position = [x1 + x2 + 1, y1 + y2 + 1] / 2; %center point
    %state.eK = eK;
    state.Q = Q;
    state.X = X;
    state.Y = Y;
    state.bins = bins;
    location = [x1, y1, state.size];
    state.B = B;
    
    %KALMAN
    k_q = ((w*h)/(width*height))*50;
      
    %matrike za NCA model gibanja
%     k_A = [1 0 1 0 (1/2) 0; 0 1 0 1 0 (1/2); 0 0 1 0 1 0; 0 0 0 1 0 1; 0 0 0 0 1 0; 0 0 0 0 0 1];
%     k_Q = k_q*[(1/20) 0 (1/8) 0 (1/6) 0; 0 (1/20) 0 (1/8) 0 (1/6); (1/8) 0 (1/3) 0 (1/2) 0;  0 (1/8) 0 (1/3) 0 (1/2); (1/6) 0 (1/2) 0 1 0; 0 (1/6) 0 (1/2) 0 1];
%     k_C = [1 0 0 0 0 0; 0 1 0 0 0 0];
    
    %matrike za NCV model gibanja
%     k_A = [1 0 1 0;0 1 0 1; 0 0 1 0; 0 0 0 1]; %T = 1
%     k_Q = k_q*[(1/3) 0 (1/2) 0; 0 (1/3) 0 (1/2); (1/2) 0 1 0; 0 (1/2) 0 1];
%     k_C = [1 0 0 0; 0 1 0 0]; % za xyx.y.

    %matrike za RW model gibanja
    k_A = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    k_Q = k_q*[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0];
    k_C = [1 0 0 0; 0 1 0 0];
    
    k_R = [1 0; 0 1];
    
    k_state = zeros(size(k_A,1),1);
    k_state(1) = state.position(1);
    k_state(2) = state.position(2);
    k_covariance = eye(size(k_A,1)).*10;
    
    state.k_A = k_A;
    state.k_Q = k_Q;
    state.k_R = k_R;
    state.k_C = k_C;
    state.k_state = k_state;
    state.k_covariance = k_covariance;
end