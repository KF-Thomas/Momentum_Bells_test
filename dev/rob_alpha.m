function alpha = rob_alpha(lambda)

alpha = exp(-2.*lambda.^2)-1+sqrt(2*pi).*lambda.*erf(sqrt(2).*lambda);

end