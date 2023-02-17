function drawElipse(x0, y0, w, h, col)
    t=-pi:0.01:pi;
    x=x0+w*cos(t);
    y=y0+h*sin(t);
    plot(x,y, col)
end

