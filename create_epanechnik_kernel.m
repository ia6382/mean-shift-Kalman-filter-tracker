function kernel = create_epanechnik_kernel(w, h, sigma)
% INPUT:
% w - width
% h - height
% sigma: sigma parameter (scalar)
% OUTPUT:
% kernel: 2-D matrix

w2 = floor(w / 2);
h2 = floor(h / 2);

[X, Y] = meshgrid(-w2:w2, -h2:h2);

X = X / max(X(:));
Y = Y / max(Y(:));

kernel = (1 - ((X/sigma).^2 + (Y/sigma).^2));
kernel = kernel / max(kernel(:));
kernel(kernel<0) = 0;

end  % endfunction