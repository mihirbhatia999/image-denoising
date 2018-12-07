function [LX] = laplace(X)
    kernel = -1 * ones(3);
    kernel(2,2) = 8; 
    LX = conv2(X, kernel, 'same');  
end
