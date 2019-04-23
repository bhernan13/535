function gradMatrix = findGradient(rawdata)

% Code Author: Madie Melcer
% Code Created: 4/20/19

% This function finds the gradient of the data given across the space
% dimension. This assumes that the vertical dimension of the matrix
% provided is space, and the horizontal dimension is time.

[m,n] = size(rawdata);

for i = 1:m-1
    for j = 1:n
        gradMatrix(i,j) = abs(rawdata(i+1,j) - rawdata(i,j));
    end
end

end