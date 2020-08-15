phi = pi/2;
Tmax = 15;
lambda = 0.8;
cycles = 1E7;

[p,event,cap] = problist(lambda,phi,Tmax);

m = zeros(1,cycles); ys = m; xs = m; ws = m; zs = m;

for runs = 1:cycles
    n = rand;
    i = 1;
    while true
        if n<p(i+1)
            m(runs) = i;
            ys(runs) = event(i,1);
            xs(runs) = event(i,2);
            ws(runs) = event(i,3);
            zs(runs) = event(i,4);
            break
        end
        i=i+1;
    end
end

% Plotting Histograms
edges = 0:cap;

subplot(2,1,1)
histogram(m,edges,'normalization','probability')
xlabel("Iterated States")
% States being each tested combination of YXWZ, beating pattern due to the
% way states are generated.

subplot(2,4,5)
histogram(ys,Tmax,'normalization','probability')
xlabel("Particles at Y")

subplot(2,4,6)
histogram(xs,Tmax,'normalization','probability')
xlabel("Particles at X")

subplot(2,4,7)
histogram(ws,Tmax,'normalization','probability')
xlabel("Particles at W")

subplot(2,4,8)
histogram(zs,Tmax,'normalization','probability')
xlabel("Particles at Z")