%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reference codes
% Authors: Xinyu Zhao
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = kawahara(v)
N = length(v);
c = v(1);
A = 0.001;
alpha = 2;
f = [A, v(2:end)];
x = [flip(f), 0, f];
w = conv(x, x);
F = w(2*N+2:2*N+2+N-1)*0.5;
index = 1:1:N;
F = c*f+0.5*F-alpha*index.^2.*f+index.^4.*f;
end