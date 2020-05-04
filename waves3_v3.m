w=[-10e6:10e6];
f0 = 10e9;
T = 1e-9;
t = [-1*T:1e-12:1*T];
pulse = zeros(size(t));
pulse(abs(t)<=T/2)=1;


a=-.5*10e-9;
b=.5*10e-9;
mod = sin(2*pi*f0.*t);

Vgt = pulse.*mod;
figure;
plot(t,Vgt);
title('Vg as a function of Time');
xlabel('t[s]');ylabel('Vg[V]');

Fs=100;
n = 2^nextpow2(length(t));
Vgw = fft(Vgt,n);
f = Fs*(0:(n/2))/n;
P = abs(Vgw/n);
figure(10);
plot(f,P(1:n/2+1));
title('Vg as a function of Angular Frequency');
xlabel('w[rad/s]');ylabel('Vg[V]');


A = 20e-3;
fc=7.5e9;
wc=2*pi*fc;
c = 3e8;
ww = wc/c:(length(t)*wc/c);
k = sqrt((ww.^2)-(pi/A)^2);
figure;
plot(ww,real(k));
title('K as a function of  angular Frequency');
xlabel('w[rad/s]');ylabel('Vg[V]');
xlim([wc/c 5*wc/c]);

vp = c/sqrt(1-(fc/f0)^2);
vg = c*sqrt(1-(fc/f0)^2);


W = 2*pi*f0;
VP = W/(sqrt(((W/c)^2)-((pi/A)^2)));

syms w;
VG = sqrt(((w/c)^2)-((pi/A)^2));
ans = diff(VG,w);
a1 = double(subs(ans,w,(2*pi*f0)));
sol = 1/a1;


z0 = vp*vg/((vp-vg)*8*f0);
z= [z0, 2*z0, 3*z0, 4*z0, 1, 200];


%for b=1:length(z)
    %modulation = sin(2*pi*f0.*(t-z(b)/vg));
    %pul = ones(size(t));
    %pul(t>=(T/2)+(z(b)/vp))=0;
    %pul(t<=(-T/2)+(z(b)/vp))=0;
%     if z(b)==1
%         t = tlim;
%     else
%     end
    %subplot(2,3,b);
    %plot(t,pul.*modulation);
    %xlim([-1e-9+(z(b)/vp) 1e-9+(z(b)/vp)]);
%end
hold off;


% tau = (-vg/100000:1:vg/100000);


hold on;
figure(8)
for b=1:length(z)
    tau=t-z(b)/vg;
    modulation = sin(2*pi*f0.*(t-z(b)/vp));
    pul = ones(size(t));
    pul(tau>=(T/2))=0;
    pul(t<=(-T/2))=0;
    Vt = pul.*modulation;
    Vw = fft(Vt);
    n=min(length(ww),length(Vw));
    k=k(1:n);
    ww=ww(1:n);
    VGw = exp(-1i*(k-(ww/vg))*z(b)).*Vw;
    VG_T = ifft(VGw);
    subplot(2,3,b);
    hold on;
    plot(tau,real(VG_T));
    title(['Vg(tau,z=',num2str(z(b)),')']);
    xlabel('tau[s]');ylabel('Vg[V]');
    hold off;

   

end
hold off;

