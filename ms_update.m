function [state, location] = ms_update(state, I, varargin)
    %colour spaces
    %I = rgb2ycbcr(I);
    %I = rgb2hsv(I);
    %I = rgb2lab(I);
    
    %parameters
    alfaMax = 0.005;
    sigma = 1;
    e = 1e-4; %eps1e-6
    cl = 1; %convergance limit = pixel
    nIterMax = 10;
    scaleFactor = 0.1;
    gama = 0.3;
    
    [height, width, ~] = size(I);
    w = state.size(1);
    h = state.size(2);
    
    %scaling - pick a size that is most similar with template
    sizeDiff = scaleFactor*state.size;
    
    same = getSimil(state, I, state.size);
    smaller = getSimil(state, I, state.size - sizeDiff);
    bigger = getSimil(state, I, state.size + sizeDiff);
    
    if ((smaller > same) && (smaller > bigger))
        w = state.size(1) - sizeDiff(1);
        h = state.size(2) - sizeDiff(2);
    end
    if ((bigger > same) && (bigger > smaller))
        w = state.size(1) + sizeDiff(1);
        h = state.size(2) + sizeDiff(2);
    end
    
    %update size - prevent oversensitive scale adaptation
    w = gama*w + (1-gama)*state.size(1);
    h = gama*h + (1-gama)*state.size(2);
    
    %get epanechnik kernel
    eK = create_epanechnik_kernel(w, h, sigma);
    wP = size(eK,2);
    hP = size(eK,1);  
    
    %coordinates for ms
    [X, Y] = meshgrid(((-wP+1)/2):((wP-1)/2), ((-hP+1)/2):((hP-1)/2)); 

    center = state.position;
    m = 2;
    iter = 0;
    while(m > cl && iter < nIterMax)
        %put eK pixels outside of image to 0
        eK = outOfBorderToZero(eK, center, wP, hP, width, height);      
        
        %extract patch
        region = get_patch(I, center, 1, [wP, hP]);
        
        %get weighted histogram from search region
        P = extract_histogram(region, state.bins, eK);

        %calculate weight for each colour bin from template and current region
        P = P/sum(P(:));
        PW = state.B.*P; %background histogram weighing
        QW = state.B.*state.Q;
        PW = PW/sum(PW(:));
        QW = QW/sum(QW(:));
        V = sqrt(QW./(PW+e));
        V = V/sum(V(:));

        %backproject into the search region
        W = backproject_histogram(region, V);
        
        %put W pixels outside of image to 0
        W = outOfBorderToZero(W, center, wP, hP, width, height); 

        %calculate mean shift vector
        xMove = sum(sum(double(X).*W))./ sum(W(:));
        yMove = sum(sum(double(Y).*W))./ sum(W(:));
        move = [xMove, yMove];

        %calculate next mean (center)
        centerNew = center + move;
        m = norm(move); %vector size

        %new iteration
        center = centerNew;
        iter = iter + 1;
        
    end
    %get region around last center (same as in the loop)
    eK = outOfBorderToZero(eK, center, wP, hP, width, height);
    region = get_patch(I, center, 1, [wP, hP]);
    P = extract_histogram(region, state.bins, eK);
    P = P/sum(P(:));
    PW = state.B.*P;
    PW = PW/sum(PW(:));
    QW = state.B.*state.Q;
    QW = QW/sum(QW(:));
    
    %Bathacharaya measure for similarity
    simil = sqrt(QW.*PW);
    simil = sum(simil(:));
    omega = ((wP*hP)/(width*height))*2000;
    k_r = omega*(1-simil);
    
    %background hist extraction
    bckgborder = 10;
    bckg = get_patch(I, center, 1, [wP, hP]+2*bckgborder);
    [bh, bw, ~] = size(bckg);
    bckgWeights = ones(bh, bw);
    bckgWeights = outOfBorderToZero(bckgWeights, center, bw, bh, width, height); 
    bckgWeights(bckgborder:bh-bckgborder, bckgborder:bw-bckgborder, :) = 0;
    B = extract_histogram(bckg, state.bins, bckgWeights);
    %B = B/sum(B(:));
    %weigh so the smallest non zero value is the largest weight
    m=min(B(B>0));
    B = min((m ./ B), 1);
    B = B/sum(B(:));
    
    %update template Q
    Qnew = P;
    Qold = state.Q;
    alfa = alfaMax*(1-simil);
    state.Q = (1 - alfa)*Qold + alfa*Qnew;
    
    %update background histogram B
    beta = alfa*10;
    state.B = (1 - beta)*state.B + beta*B;
    
    %KALMAN
    [k_state,k_covariance] = kalman_update(state.k_A, state.k_C, state.k_Q, k_r*state.k_R, center', state.k_state, state.k_covariance);
    sCenterx = k_state(1);
    sCentery = k_state(2);
    
    %update struct  
    state.size = [w h];
    state.k_state = k_state;
    state.k_covariance = k_covariance;
    state.position = [sCenterx, sCentery];
%     state.position = centerNew;
    location = [state.position - state.size / 2, state.size];
end
