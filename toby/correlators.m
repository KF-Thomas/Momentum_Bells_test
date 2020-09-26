function [gYW,gXZ,gYZ,gXW] = correlators(lambda,phi,Tmax,cycles,quantum_efficiency,dark_rate,blur,R,phys,Y_bias,X_bias)
    colist_test = coordinatelist(lambda,phi,Tmax,cycles,quantum_efficiency,dark_rate,blur,R,phys,Y_bias,X_bias);
    total_pairs = 0;
    FYW = 0; FXZ = 0; FYZ = 0; FXW = 0;
    FYW2 = zeros(1,cycles); FXZ2 = zeros(1,cycles); FYZ2 = zeros(1,cycles); FXW2 = zeros(1,cycles);
    for i = 1:cycles
        n = size(colist_test{i},1);
        no_pairs = (n^2-n)/2;
        total_pairs = total_pairs + no_pairs;
        if size(colist_test{i})==0
            continue
        end
        no_Ys = length(colist_test{i}(colist_test{i}(:,2)>0 & colist_test{i}(:,3)>0));
        no_Xs = length(colist_test{i}(colist_test{i}(:,2)>0 & colist_test{i}(:,3)<0));
        no_Ws = length(colist_test{i}(colist_test{i}(:,2)<0 & colist_test{i}(:,3)>0));
        no_Zs = length(colist_test{i}(colist_test{i}(:,2)<0 & colist_test{i}(:,3)<0));
    
        FYW = FYW + no_Ys * no_Ws;
        FXZ = FXZ + no_Xs * no_Zs;
        FYZ = FYZ + no_Ys * no_Zs;
        FXW = FXW + no_Xs * no_Ws;
        
        FYW2(i) = no_Ys * no_Ws / no_pairs;
        FXZ2(i) = no_Xs * no_Zs / no_pairs;
        FYZ2(i) = no_Ys * no_Zs / no_pairs;
        FXW2(i) = no_Xs * no_Ws / no_pairs;
    end
    
    fYW = FYW / total_pairs;
    fXZ = FXZ / total_pairs;
    fYZ = FYZ / total_pairs;
    fXW = FXW / total_pairs;
    
    fYW2 = mean(FYW2(find(FYW2)));
    fXZ2 = mean(FXZ2(find(FXZ2)));
    fYZ2 = mean(FYZ2(find(FYZ2)));
    fXW2 = mean(FXW2(find(FXW2)));
    
    colist2 = [];
    for i=1:cycles
        colist2 = [colist2; colist_test{i}];
    end
    n = size(colist2,1);
    no_pairs = (n^2-n)/2;
    no_Ys = length(colist2(colist2(:,2)>0 & colist2(:,3)>0));
    no_Xs = length(colist2(colist2(:,2)>0 & colist2(:,3)<0));
    no_Ws = length(colist2(colist2(:,2)<0 & colist2(:,3)>0));
    no_Zs = length(colist2(colist2(:,2)<0 & colist2(:,3)<0));

    SYW = no_Ys * no_Ws;
    SXZ = no_Xs * no_Zs;
    SYZ = no_Ys * no_Zs;
    SXW = no_Xs * no_Ws;
    
    sYW = SYW / no_pairs;
    sXZ = SXZ / no_pairs;
    sYZ = SYZ / no_pairs;
    sXW = SXW / no_pairs;

    gYW = fYW/sYW;
    gXZ = fXZ/sXZ;
    gYZ = fYZ/sYZ;
    gXW = fXW/sXW;
end