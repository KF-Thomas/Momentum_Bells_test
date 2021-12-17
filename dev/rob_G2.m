function g2 = rob_G2(lambda,h,theta)
g2 = h(3) + h(3).*h(1)./16.*(1+cos(theta)).*prod(rob_alpha(lambda.*h(2))).*prod((lambda.*h(2)).^(-2));
end