%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reference codes
% Authors: Xinyu Zhao
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% solve for zero
A = 0.001/2;
alpha = 2;
c = 1;
N = 100;
v0 = [c,zeros(1,N-1)];
%v0(2) = 0.001;
v = fsolve(@kawahara, v0);

fHat = [0,A,v(2:end),flip(v(2:end))];
f = real(ifft(fHat)*N*2);
plot(0:2*pi/(2*N):2*pi-2*pi/(2*N), f);


% verify
c = v(1);
index = [0:N, -N+1:-1];
ddfHat = -index.^2.*fHat;
ddddfHat = -index.^2.*ddfHat;
FHat = c*fHat+0.5*fHat.*fHat+alpha*ddfHat+ddddfHat;
FHat(1) = 0;
F = real(ifft(FHat)*(2*N));
norm(F,1)/(2*N);

v0 = v;
error = [];
figure
for A = 0.001:0.001:0.01
    alpha = 2;
    N = 100;
    %v(2) = A;
    fHat = [0,A,v(2:end),flip(v(2:end))];
    f = real(ifft(fHat)*N*2);
    plot(0:2*pi/(2*N):2*pi-2*pi/(2*N), f);
    hold on

    c = v(1);
    index = [0:N, -N+1:-1];
    ddfHat = -index.^2.*fHat;
    ddddfHat = -index.^2.*ddfHat;
    FHat = c*fHat+0.5*fHat.*fHat+alpha*ddfHat+ddddfHat;
    FHat(1) = 0;
    F = real(ifft(FHat)*100);
    error = [error, norm(F,1)/(2*N)];
    v0 = v;
end

hold off

figure
semilogy(0.001:0.001:0.01, error, '*-');
hold off