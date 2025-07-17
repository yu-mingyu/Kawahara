%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Codes for numerical solutions to the Kawahara equation
% Authors: Mingyu Yu, Linjie Ying, Juanita Gasca, Chiara L. Carnevale
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
alpha=5;
c=4;
f0 =0;
f1 = 0.001/2;

% f(x) = 0.001cos(x) as initial guess 

% assume a is 2 sided, from a_{-N} to a_{N}
% Assume N=2^m-1. 
N=2^6-1;
% This will be helpful when using verifyfft, which requires things to have
% a length which is a power of 2.
a = zeros(N-1,1);
a(1)=f1;
% 
for j=1:50
    F_local = F(a,c,alpha,f0,f1);

    DF_local = DF(a,c,alpha,f0,f1);

    dT = DF_local\F_local;
    
    if any(isnan (dT))   || norm(dT) > 1e5
        % If something goes wrong, we break the loop
        break
    elseif norm(dT) <1e-14
        % If the algorithm seems to have converged, we break the loop
        break
    end
    a = a-dT(2:end,:); % update a
    c= c-dT(1); % update c
end

% F_local = F(intval(a),alpha,beta,gamma,epsilon);

Plot_Fourier(a,f0,f1);
a
c
Plotf(a,f0,f1, pi);
norm(F(a,c,alpha,f0,f1),1)
return




function DF_out = DF(a,c,alpha,f0,f1)
    
    %input a must go from f2 to fN
    N = length(a)+1;
    
    %we create the right lower diagonal submatrix of M1
    K = (2:N)';
    diag_M = diag(c-alpha.* K.^2 + K.^4);

    %we add the first column and row of M1 manually
    M1 = [a,diag_M];
    M1 = [zeros(1,N);M1];
    M1(1,1) = f1;
    
    %we create M2
    col1 = [f0;f1;a(1:end-1)];
    row1= col1';
    M2 = toeplitz(col1,row1);
    M2 = [zeros(N,1),M2(:,2:end)];

    DF_out = M1 + M2;
end

function out = F(a,c,alpha,f0,f1)

    % assume a is from a_{2} to a_{N}
    N = (length(a)+1); % N=2^m-1;
    
    % Define the Derivative operator in Fourier
    K = (2:N)';
    
    % We compute the principal linear part
    Lin = c.*a - alpha.*(K.^2).*a + (K.^4).*a;
    
    % We compute the nonlinearity 
    [spec_out ] = convolution_product(a,2,f0,f1);

    %Here we truncate the convolution to get from f2 to fN
    conv_trunc = spec_out(N+3:end);

    % Sum everything together
    out = Lin + (1/2) .*conv_trunc;

    %add the first component of F manually
    F1 = c*f1 + (1/2).* spec_out(N+2) - alpha*f1 + f1;
    out = [F1;out];
end

function [spec_out ] = convolution_product(a,power,f0,f1)
    %for fft we need the vector to go from f(-N) to f(N)
    % input a must be f_2 ... f_N

    a_flip = flipud(a);
    conv_a = [a_flip ; f1 ; f0 ; f1 ; a];

    % Assume the indices of conv_a are ordered like: 
    %  -2    -1     0     1     2
    % And that a_k corresponds to the coefficient of e^{ikx}
    N = (length(conv_a)-1)/2; % N=2^m-1;
    n = (length(conv_a)+1)/2;

    % Adds in an extra zero corresponding to a_{-(N+1)} 
    % We ideally want the length of a_pad to be a power of 2.
    a_pad = [zeros(3*n,1);0;conv_a;zeros(3*n,1)];
    
    % Rearrange so that the coefficients are ordered like:
    % 0     1     2    -3    -2    -1
    a_shift = fftshift(a_pad);

    % Normalize things for ifft
    a_shift_normalized = a_shift*length(a_shift);

    % Go to grid space
    a_grid = ifft(a_shift_normalized);
    
    % Apply nonlinearity
    a_nonlinear = a_grid.^power;
    
    % Return to coefficient space. 
    spec_all = fft(a_nonlinear );

    % Renormalize
    spec_out = spec_all/length(spec_all);

    % Go back to coefficient ordering of :
    % -3    -2    -1     0     1     2
    spec_out  = fftshift(spec_out);
    % Truncate, so that we run from a_{-N} to a_{N} .
    spec_out = spec_out(3*n+1:3*n+2*n);
    spec_out = spec_out(2:end);

end


% I don't know what does this plot mean Ask TA lol

function a_grid = Plot_Fourier(a,f0,f1)
    % Adds in an extra zero corresponding to a_{-(N+1)} 
    a_flip = flipud(a);
    conv_a = [a_flip ; f1 ; f0 ; f1 ; a];
    a_pad = [0;conv_a];

    % Rearrange so that the coefficients are ordered like:
    % 0     1     2    -3    -2    -1
    a_shift = fftshift(a_pad);

    % Normalize things for ifft
    a_shift_normalized = a_shift*length(a_shift);

    % Go to grid space
    a_grid = ifft(a_shift_normalized); 

    plot(real(a_grid))
end


function plotf = Plotf(a, f0, f1, M)
    x = linspace(-M,M,M*100);
    N = length(a);
    y = 2 * f1 * cos(x);
    for i = 1:N
        y = y + 2*a(i)*cos((i+1)*x);
    end
    plot(x,y);
end
