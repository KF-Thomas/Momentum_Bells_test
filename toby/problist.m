% Given a mode occupancy, an interferometer phase and a cap on the number
% of detections, the following function will divide [0 1] into a section
% for each potential state, with the size determined by the probability of
% that state.  Generating a random number will sit inside one of these
% sections and hence correspond to a specific detection state.
%% Function:
function [p,event] = problist(lambda,phi,Tmax,phys,Y_bias,X_bias)
    cap = 1+13*(Tmax)/6 + (3/2)*Tmax^2 + Tmax^3/3; % The number of possible states given 
    problist = zeros(1,cap+1);
    YXWZ = zeros(cap,4);
    index = 1;
    for A = 0:Tmax
        for C = 0:Tmax
            for B = max(C-A,0):(Tmax-A)
                index = index+1;
                problist(index) = problist(index-1)+cyxwz(A,B,C,A+B-C,lambda,phi,phys,Y_bias,X_bias);
                YXWZ(index-1,:) = [A B C A+B-C];
            end
        end
    end
    problist(end) = 1;
    p = problist;
    event = YXWZ;
    cap = cap;
end