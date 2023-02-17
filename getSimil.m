function simil = getSimil(state, I, dim)
    %parameters
    sigma = 1;
    e = 1e-4; %eps
    cl = 1; %convergance limit = pixel
    nIterMax = 10;

    %get dimensions of image and region
    [height, width, ~] = size(I);
    w = dim(1);
    h = dim(2);
    if(w < 1 || h < 1)
        simil = -1;
    	return;
    end
    
    %get epanechnik kernel
    eK = create_epanechnik_kernel(w, h, sigma);
    w = size(eK,2);
    h = size(eK,1);  
    
    %coordinates for ms
    [X, Y] = meshgrid(((-w+1)/2):((w-1)/2), ((-h+1)/2):((h-1)/2)); 

    center = state.position;
    m = 2;
    iter = 0;
    while(m > cl && iter < nIterMax)
        %put eK pixels outside of image to 0
        eK = outOfBorderToZero(eK, center, w, h, width, height);      
        
        %extract patch
        region = get_patch(I, center, 1, [w, h]);
        
        %get weighted histogram from search region
        P = extract_histogram(region, state.bins, eK);

        %calculate weight for each colour bin from template and current region
        P = P/sum(P(:));
        PW = state.B.*P; %background histogram weighing
        PW = PW/sum(PW(:));
        QW = state.B.*state.Q;
        QW = QW/sum(QW(:));
        V = sqrt(QW./(PW+e));
        V = V/sum(V(:));

        %backproject into the search region
        W = backproject_histogram(region, V);
        
        %put W pixels outside of image to 0
        W = outOfBorderToZero(W, center, w, h, width, height); 

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
    eK = outOfBorderToZero(eK, center, w, h, width, height);
    region = get_patch(I, center, 1, [w, h]);
    P = extract_histogram(region, state.bins, eK);
    P = P/sum(P(:));
    PW = state.B.*P;
    PW = PW/sum(PW(:));
    QW = state.B.*state.Q;
    QW = QW/sum(QW(:));

    %calculate Batacharaya measure
    simil = sqrt(QW.*PW);
    simil = sum(simil(:));
end