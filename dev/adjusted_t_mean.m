function tc = adjusted_t_mean(t0,t1)
global const
t0 = t0+19.5-20.854;
t1 = t1+19.5-20.854;
d = const.fall_distance;
g0 = const.g0;
tf = sqrt(2*d/g0);%arrival time of zero velocity particles
tavg = (t0+t1)./2;
A = 1/2.*tavg - tf.^2 .*tavg./(t0.*t1);
tc = A + sqrt(A.^2+2.*tf.^2);
tc = tc - (19.5-20.854);
end