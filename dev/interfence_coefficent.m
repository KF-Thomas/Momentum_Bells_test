function val = interfence_coefficent(Y,X,W,Z,phi,lmd)
%gives the coefficents of th state after interference for a give phi and
%lmd

val = 0 ;
if Y+W == X+Z
for m = 0:(W+Y)
    for l = max([m-Y,0]):min([m,W])
        for n = max([m-X,0]):min([m,Z])    
            val = val + (1-lmd^2)*(lmd/2)^(W+Y)*sqrt(factorial(Y)*factorial(X)*factorial(W)*...
                factorial(Z))*1/(factorial(m).*factorial(W+Y-m)).*1i^(W+Z+2*m-2*l-2*n)*...
                exp(-2*1i*phi*m)*nchoosek(W+Y-m,W-l)*nchoosek(W+Y-m,Z-n)*nchoosek(m,l)*...
                nchoosek(m,n);
        end
    end
end
end
end