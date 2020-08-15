function b = bincof(n,k)
    b = factorial(abs(n))./(factorial(abs(n-k)).*factorial(abs(k)));
end
    