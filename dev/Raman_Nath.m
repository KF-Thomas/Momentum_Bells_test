function dydt = Raman_Nath(t,c,N,wr,b,k)
%VDP1  Evaluate the van der Pol ODEs for mu = 1
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

modes = 2*N-1;
ac = 7.933e-15;
k0 = sqrt(wr/ac);
% wr = 1;
% Omega = 1;

dydt = zeros(modes,1);

for ii = 2:(modes-1)
    n = ii-N;
    dydt(ii,1) = -1i*(ac*(2*n*k0+k)^2*c(ii)+Omega(t,b)/2*(exp(-1i*theta(t,b))*c(ii-1)+exp(1i*theta(t,b))*c(ii+1)));
end
if N>1
    n = -(N-1);
    dydt(1,1) = -1i*(ac*(2*n*k0+k)^2*c(1)+Omega(t,b)/2*(exp(1i*theta(t,b))*c(2)));
    n = N-1;
    dydt(end,1) = -1i*(ac*(2*n*k0+k)^2*c(end)+Omega(t,b)/2*(exp(-1i*theta(t,b))*c(end-1)));
end

end

function amp = Omega(t,b)
% amp = gaussian_pulse(t,b(1),b(2)).*b(3);
% amp = super_gaussian_pulse(t,b(1),b(2),b(4)).*b(3);
amp = sinc_pulse(t,b(1),b(2)).*b(3);
% amp = square_pulse(t,b(1),b(2)).*b(3);
% amp = sinc_pulse(t,b(1),b(2)).*gaussian_pulse(t,b(1),b(6)).*b(3);
% amp = sinc_pulse(t,b(1),b(2)).*b(3);%.*cos(pi*(t-b(1))/(2*b(1))).^2
% amp = sinc_pulse(t,b(1),b(2)).*b(3).*cos(pi*(t-b(1))/(2*b(1))).^2;
% amp = hamming_pulse(t,b(1),b(2)).*b(3);%.*cos(pi*(t-b(1))/(2*b(1))).^2
% amp = sqrt(abs(sinc_pulse(t,b(1),b(2)).*b(3).*cos(pi*(t-b(1))/(2*b(1))).^2));
end

function phase = theta(t,b)
% phase = 0;
phase = (b(4)*t+b(5));%b(4)*t;
% phase = (b(5)*t+b(6)*t^2);%b(4)*t;
end

function amp = square_pulse(t,start_time,duration)
if t>start_time && t<start_time+duration
    amp = 1;
else
    amp = 0;
end
end

function amp = gaussian_pulse(t,t0,alpha)
amp = exp(-alpha*(t-t0)^2);
end

function amp = super_gaussian_pulse(t,t0,alpha,n)
amp = exp(-alpha*(t-t0).^n);
end

function amp = sinc_pulse(t,t0,alpha)
amp = sinc((t-t0)/alpha);
end

function amp = hamming_pulse(t,t0,L)
if abs(t-t0)>L/2
    amp = 0;
else
    amp = cos(pi*(t-t0-L/2)/(L)).^2;
end
end