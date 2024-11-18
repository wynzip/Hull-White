function paramsFFT = parametersFFT(M, dz)

% Setting up the parameters of the FFT grid using the known constraints,
% starting from two fixed parameters (M,dz)
%
% INPUT
% M         integer used to obtain N
% dz        step of the integration grid
%
% OUTPUT
% paramsFFT struct containing all the elements below:
% N         number of points in each grid
% z1        lowest point of integration grid (truncation)
% zGrid     integration grid (truncation)
% x1        lowest point of moneyness grid
% xGrid     moneyness grid
% dz        step of the integration grid
% dx        step of the moneyness grid
%

    % Evaluate number of points in each grid
    N = 2^M; 
    
    % Truncation of the integral grid (integration variable csi)
    z1 = -dz*(N-1)/2;
    zN = -z1;
    zGrid = (z1:dz:zN)';
    
    % Moneyness grid for computation
    dx = 2*pi/(dz*N); % known constraint: dx*dz = 2*pi/N
    x1 = -dx*(N-1)/2;
    xN = -x1;
    xGrid = (x1:dx:xN)';
    
    % Storing inside the struct all the parameters found
    paramsFFT.N = N;
    paramsFFT.z1 = z1;
    paramsFFT.zGrid = zGrid;
    paramsFFT.x1 = x1;
    paramsFFT.xGrid = xGrid;
    paramsFFT.dz = dz;
    paramsFFT.dx = dx;

end % parametersFFT