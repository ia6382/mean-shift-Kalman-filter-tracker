function A = outOfBorderToZero(A, center, w, h, width, height)
    %w, h = dimensions of pathc
    %widht, height dimesnions of image
    
    %get left upper corner and lower right from center point and size
	x1 = floor(center(1) - (w-1)/2);
    y1 = floor(center(2) - (h-1)/2);
    x2 = floor(center(1) + (w-1)/2);
    y2 = floor(center(2) + (h-1)/2);
    %calculate how many pixels fell out of image
    if(x1 < 1)
        out = -x1+1;
        A(:, 1:1+out) = 0;
    end
    if(x2 > width)
        out = width-x2;
        A(:, w+out:w) = 0;
    end
    if(y1 < 1)
        out = -y1+1;
        A(1:1+out, :) = 0;
    end
    if(y2 > height)
        out = height-y2;
        A(h+out:h, :) = 0;
    end
end

