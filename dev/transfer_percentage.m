function trans_percentage = transfer_percentage(bs,ord,opts,k)
trans_percentage=[];
if opts.control == 1
    for ii = 1:length(bs)
        b=bs(ii,:);
        [t,y] = Raman_Nath_Solver(b,opts,k);
        trans_percentage(ii,:) = y(end,ord);
    end
else
    trans_percentage=zeros(length(k),length(ord));
    for ii = 1:length(k)
        k_current = k(ii);
        [t,y] = Raman_Nath_Solver(bs,opts,k_current);
        trans_percentage(ii,:) = y(end,ord);
    end
end
end