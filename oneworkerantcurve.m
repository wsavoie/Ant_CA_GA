figure(25)
hold on;
for L=[1]
a = 0.001:.001:.3; % alpha
g = 1/25; % "switching rate" of excavation 
v = .75; % Speed of ant (1 cell per frame)
L = 6; % Length of tunnel
S = .9; % Switching rate due to reversal
p = 2*((L/2-.5*(1./a-v/g))./(.5*(1./a-v/S))+2)/v.*(1./(v*a)+2*L/v+1/S.*((L-.5.*(1./a-v/g))./(.5*(1./a-v/S))+1)+1/g).^-1;
% p is rho, or density of tunnel

q = (1./(v*a)+2*L/v+1/S.*((L-.5.*(1./a-v/g))./(.5*(1./a-v/S))+1)+1/g).^-1;
% q is flow rate, in ants per frame

conv = 0.5; % conversion, 0.5 seconds per frame.

plot(p,q/conv)
xlabel('\rho');
ylabel('q');
end