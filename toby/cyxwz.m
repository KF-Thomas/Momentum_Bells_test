function c = cyxwz(y,x,w,z,T,lambda,phim,phip)
    alpha = (1-lambda.^2).*exp(1i*(phip*z+phim*w)).*(lambda./2).^T*sqrt(factorial(y)*factorial(z)*factorial(w)*factorial(x));
    psi = 0;
    for k = 0:T
        %runk = k
        chi = 0;
        tkfact = factorial(T-k);
        for l = max(k-y,0):min(k,w)
            %runl = l
            epsilon = 0;
            wlfact = factorial(w-l);
            for n = max(k-x,0):min(k,z)
                %runn = n
                znfact = factorial(z-n);
                beta = 1i.^(w+z+2.*k-2.*l-2.*n);
                gamma = exp(-1i.*(phip+phim).*k);
                delta = (tkfact./wlfact).*(1./factorial(abs(T-k-w+l))).*(tkfact./znfact).*(1./factorial(abs(T-k-z+n))).*nchoosek(k,l).*nchoosek(k,n);
                epsilon = epsilon + beta.*gamma.*delta;
            end
            chi = chi+epsilon;
        end
        psi = psi+(factorial(T-k).*factorial(k)).^(-1).*chi;
    end
    c = abs(alpha.*psi).^2;
end
