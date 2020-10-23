function Y = mypdist(X)
    [n,~] = size(X);
    
    if n < 2
        Y = zeros(1,0);
        return;
    end
    
    north = X(X(:,3)>0,:);
    south = X(X(:,3)<0,:);
    northeast = north(north(:,2)>0,:);
    northwest = north(north(:,2)<0,:);
    
    [Nn,~] = size(north);
    [Sn,~] = size(south);
    [NEn,~] = size(northeast);
    [NWn,~] = size(northwest);
    
    Y = zeros(1,Nn*Sn);
    k = 1;
    
    
    for i = 1:NEn
        dsq = (-northeast(i,1) - south(1:Sn,1)).^2 + (-northeast(i,2) - south(1:Sn,2)).^2 + (-northeast(i,3)-south(1:Sn,3)).^2;
        Y(k:(k+Sn-1)) = sqrt(dsq);
        k = k + Sn;
    end
    for i = 1:NWn
        dsq = (-northwest(i,1) - south(1:Sn,1)).^2 + (-northwest(i,2) - south(1:Sn,2)).^2 + (-northwest(i,3)-south(1:Sn,3)).^2;
        Y(k:(k+Sn-1)) = -sqrt(dsq);
        k = k + Sn;
    end
    for i = 1:Nn-1
        dsq = (-north(i,1)-north((i+1):Nn,1)).^2 +(-north(i,2)-north((i+1):Nn,2)).^2 + (2-north(i,3)-north((i+1):Nn,3)).^2;
        Y(k:(k+Nn-i-1)) = 10+sqrt(dsq);
        k = k+Nn-i;
    end
    
    for i = 1:Sn-1
        dsq = (-south(i,1)-south((i+1):Sn,1)).^2 +(-south(i,2)-south((i+1):Sn,2)).^2 + (-2-south(i,3)-south((i+1):Sn,3)).^2;
        Y(k:(k+Sn-i-1)) = -10+sqrt(dsq);
        k=k+Sn-i;
    end

end