clc
lambda = 0.3;
phi = pi/2;
Tmax = 7;
cap = (Tmax+1)/6 + (Tmax+1)^3*5/6+1;
prob = [];
picklist = zeros(1,cap);
YXWZ = zeros(cap,4);
index = 1;
runs = [0];
trial = [0];
cs = [0];
for A = 0:Tmax
    for B = 0:Tmax
        for C = max(B-A,0):Tmax
            index = index+1;
            picklist(index) = picklist(index-1)+cyxwz(A,B,C,A+C-B,A+C,lambda,0,phi);
            YXWZ(index,:) = [A B C A+C-B];
        end
    end
end


n = rand;
%n = rand(1,20);
m = zeros(1,20);
i=1;
while n<picklist(i)
    i = i+1;
end
i
return

for j=1:length(n);
    i = 1;
    while n(j) < picklist(i)
        i=i+1;
    end
    m(j) = picklist(i-1);
end
return
picklist(i-1)
picklist(i)
return
first = picklist(1)
last = picklist(end)