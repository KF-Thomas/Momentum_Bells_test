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