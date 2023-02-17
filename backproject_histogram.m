function projection = backproject_histogram(image, histogram)
% INPUT:
% image: image patch - must be colour image (RGB or different colorspaces)
% histogram: 3-dimmensional matrix (obtained with extract_histogram)
% OUTPUT:
% projection: histogram backprojection into the image patch

[h, w, c] = size(image);

if (ndims(histogram) ~= c)
    error('Dimensions do not match');
end;

bins = size(histogram) - 1;

if isa(image, 'uint8')
    image = round((reshape(double(image), [h * w, c]) ./ 255.0) .* repmat(bins, h * w, 1)) + 1;
elseif isa(image, 'double')
    image = round((reshape(max(min(image, 1), 0), [h * w, c])) .* repmat(bins, h * w, 1)) + 1;
else
    error('Unsupported image format');
end

channels = mat2cell(image, h * w, ones(c, 1));

projection = reshape(histogram(sub2ind(size(histogram), channels{:})), h, w);
