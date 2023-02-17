function histogram = extract_histogram(image, bins, weights)
% INPUT:
% image: image patch - must be colour image (RGB or different colorspaces)
% bins: number of bins (scalar)
% weights: weighting kernel (the same size as input image patch)
% OUTPUT:
% histogram: 3-dimmensional matrix


if nargin < 2
    bins = 8;
end;

if nargin < 3
    weights = [];
end;

if (size(image, 3) ~= 3)
    error('image_histogram:numberOfChannels', 'Input image must contain three channels.')
end

if numel(bins) ~= 3
    bins = [bins(1), bins(1), bins(1)];
end;

n = size(image,1) * size(image,2);

if isa(image, 'uint8')
    image = round((reshape(double(image), [n, 3]) ./ 255.0) .* repmat(bins - 1, n, 1)) + 1;
elseif isa(image, 'double')
    image = round((reshape(double(image), [n, 3])) .* repmat(bins - 1, n, 1)) + 1;
else
    error('Unsupported image format');
end

if ~isempty(weights)
    histogram = accumarray(image, weights(:), bins);
else
    histogram = accumarray(image, 1, bins);
end;