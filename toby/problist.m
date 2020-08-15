function [p,event,cap] = problist(lambda,phi,Tmax)
    cap = 2+7*(Tmax)/3 +2*Tmax^2+2*Tmax^3/3;
    problist = zeros(1,cap+1);
    YXWZ = zeros(cap,4);
    index = 1;
    for A = 0:Tmax
        for B = 0:Tmax
            for C = max(B-A,0):min(Tmax,Tmax+B-A)
                index = index+1;
                problist(index) = problist(index-1)+cyxwz(A,B,C,A+C-B,lambda,phi);
                YXWZ(index-1,:) = [A B C A+C-B];
            end
        end
    end
    problist(end) = 1;
    p = problist;
    event = YXWZ;
    cap = cap;
end