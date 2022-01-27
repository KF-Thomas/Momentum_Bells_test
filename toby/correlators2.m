function [g_side,g_diag,g,f,s] = correlators2(lambda,phi,Tmax,cycles,quantum_efficiency,dark_rate,blur,dr,phys,Y_bias,X_bias)
    colist_test = coordinatelist(lambda,phi,Tmax,cycles,quantum_efficiency,dark_rate,blur,phys,Y_bias,X_bias);
    
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
    dr2 = 1+sqrt(1-(dr^2)/4);
    edges = [0 dr dr2 2.5 3+dr2 5.5];
    %edges = 0:0.1:5.5;
    h1 = histogram(dlist,edges,'Normalization','Probability');
    f = h1.Values;
    h2 = histogram(d2list,edges,'Normalization','Probability');
    s = h2.Values;
    g = (f)./(s);
    g_side = g(end);
    g_diag = g(3);
end