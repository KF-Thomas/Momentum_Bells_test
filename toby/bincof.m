%% Function to mimic the choose function or binomial coefficients (used in the CYXWZ file)
function b = bincof(n,k)
    b = factorial(abs(n))./(factorial(abs(n-k)).*factorial(abs(k)));
end
    