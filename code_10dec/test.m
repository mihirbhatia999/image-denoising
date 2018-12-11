X = rand(5,5);
[h,w] = size(X);
delta = 5;

qwer = (power(norm((laplace(X)),'fro'),2))+(w*delta*delta);

tyu = trace(transpose(laplace(X))*laplace(X) + delta*delta*eye(w)) ;

qwer - tyu 

