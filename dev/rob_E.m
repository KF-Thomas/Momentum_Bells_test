function E_val = rob_E(lambda,h,theta)
E_val=cos(theta).*h(1).*prod(rob_alpha(lambda))./(16.*prod(lambda.^2)+h(1).*prod(rob_alpha(lambda)));
end