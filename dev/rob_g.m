function g2 = rob_g(lambda,h,theta)
g2 = 1 + h(3) + h(1)./16.*(1+cos(theta)).*prod(rob_alpha(lambda.*h(2))).*prod((lambda.*h(2)).^(-2));
end