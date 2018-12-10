%% convolution of X and kernel L
%padding around the image X to preserve size
function [LX] = laplace(X)
    fprintf('laplace !!')
    
    [h,w] = size(X);
    L = -1 * ones(3);
    size(L)
    class(L)
    L(2,2) = 8.0;
    
    image = [zeros(h,1), X , zeros(h,1)]; 
    image = [zeros(1,w+2); image; zeros(1,w+2)];
    padded_X = image ;
    new_image = X;
    
    %now doing the convolution on new image
    for i = 1:h
        for j = 1:w
            new_image(i,j) = sum(sum((padded_X(i : i +2 ,j :j+2 ).* L),2));
        end
    end
    LX = new_image;
end