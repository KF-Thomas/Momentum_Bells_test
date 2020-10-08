function [g_side,g_diag,g] = correlators2(lambda,phi,Tmax,cycles,quantum_efficiency,dark_rate,blur,R,phys,Y_bias,X_bias)
    colist_test = coordinatelist(lambda,phi,Tmax,cycles,quantum_efficiency,dark_rate,blur,R,phys,Y_bias,X_bias);
    %total_pairs = 0;
    %FYW = 0; FXZ = 0; FYZ = 0; FXW = 0;
    %FYW2 = zeros(1,cycles); FXZ2 = zeros(1,cycles); FYZ2 = zeros(1,cycles); FXW2 = zeros(1,cycles);
    dlist = [];
    for i = 1:cycles
        temp = colist_test{i};
        dlist = [dlist mypdist(temp)];
    end
    
    colist2 = [];
    for i=1:cycles
        colist2 = [colist2; colist_test{i}];
    end
    d2list = mypdist(colist2);
    edges = 0:0.01:2+3*blur;
    h1 = histogram(dlist,edges,'Normalization','Probability');
    f = h1.Values;
    h2 = histogram(d2list,edges,'Normalization','Probability');
    s = h2.Values;
    g = (f)./(s);
    g_side = g(1);
    g_diag = g(end);
end