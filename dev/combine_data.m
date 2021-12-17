function [x_out,y_out,w_out] = combine_data(x,y,w)

v = unique( x, 'stable' );

if nargin == 3
    x_out = zeros(length(v),1);
y_out = zeros(length(v),1);
w_out = zeros(length(v),1);
    for ii = 1:length(v)
        cv = v(ii);
        cv_index = (x==cv);
        x_out(ii) = cv;
        y_out(ii) = nansum(1./abs(w(cv_index)).*y(cv_index,:))./nansum(1./abs(w(cv_index)));
        w_out(ii) = nanmean(w(cv_index))./sqrt(sum(cv_index));
    end
else
    for ii = 1:length(v)
        cv = v(ii);
        cv_index = (x==cv);
        x_out(ii) = cv;
        y_out(ii,:) = mean(y(cv_index,:),1);
    end
end

end