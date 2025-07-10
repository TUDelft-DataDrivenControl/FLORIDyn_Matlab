%     dir_mode = Direction_Constant()
%     WindDir = 270
%     iT = [1, 2, 3]
%     phi = getWindDirT(dir_mode, WindDir, iT, nothing)
%     for ph in phi
%         @test ph ≈ 270.0
%     end

clear all
addpath('./WindField/Direction_Constant');
WindDir = 270;
iT = [1, 2, 3];
phi = getWindDirT(WindDir, iT);

tol = 1e-8; % set a tolerance for floating-point comparison
for i = 1:length(phi)
    ph = phi(i);
    assert(abs(ph - 270.0) < tol);
end
rmpath('./WindField/Direction_Constant')

%     dir_mode = Direction_Constant_wErrCov()
%     struct WindDirType
%         Data::Float64
%         CholSig::Matrix{Float64}
%     end
%     WindDir = WindDirType(270.0, cholesky(Matrix{Float64}(I, 3, 3)).L)
%     iT = [1, 2, 3]
%     phi = getWindDirT(dir_mode, WindDir, iT)
%     result = [269.6402710931765, 271.0872084924286, 269.5804103830612]
%     for (i, ph) in pairs(phi)
%         @test ph ≈ result[i]
%     end

clear all
addpath('./WindField/Direction_Constant_wErrorCov');
rng(1234);

% Cholesky factor of 3x3 identity
L = chol(eye(3), 'lower');

% WindDir structure
WindDir.Data = 270.0;
WindDir.CholSig = L;

iT = [1, 2, 3]';

% Expected result
result = [269.0527533560426, 270.54014974707036, 269.78339785902375];

phi = getWindDirT(WindDir, iT);
for i = 1:length(phi)
    assert(abs(phi(i) - result(i)) < 1e-8, 'Test failed at index %d', i);
end
rmpath('./WindField/Direction_Constant_wErrorCov');

%   dir_mode = Direction_EnKF_InterpTurbine()
%   # Suppose WindDir is a matrix where each row is [time, phi_T0, phi_T1, ...]
%     WindDir = [
%         0.0  10.0  20.0
%         1.0  12.0  22.0
%         2.0  14.0  24.0
%     ]
%     phi = getWindDirT_EnKF(dir_mode, WindDir, 1, 0.5)
%     @test phi ≈ 11.0

clear all
addpath('./WindField/Direction_EnKF_InterpTurbine');
WindDir = [
    0.0  10.0  20.0
    1.0  12.0  22.0
    2.0  14.0  24.0
];
phi = getWindDirT_EnKF(WindDir, 1, 0.5);
assert(abs(phi - 11.0) < 1e-8, 'Test failed');
rmpath('./WindField/Direction_EnKF_InterpTurbine');

% dir_mode = Direction_Interpolation_wErrorCov()
% 
% # Example wind direction data (time, phi)
% wind_data = [
%     0.0  10.0;
%     5.0  20.0;
%     10.0 30.0
% ]
% 
% # Example Cholesky factor (for 2 turbines)
% chol_sig = cholesky([1.0 0.5; 0.5 1.0]).L
% 
% # Create WindDir instance
% WindDir = WindDirMatrix(wind_data, chol_sig)
% 
% # Example turbine indices (for 2 turbines)
% iT = [1, 2]
% 
% # Example time
% t = 7.0
% 
% # Call the function
% phi = getWindDirT(dir_mode, WindDir, iT, t)
% @test size(phi) == (2,1)
% @test phi[1] ≈ 24.92903352636283
% @test phi[2] ≈ 24.363944731838128

clear all
addpath('./WindField/Direction_Interpolation_wErrorCov');

% Wind direction data: [time, direction]
wind_data = [
    0.0, 10.0;
    5.0, 20.0;
    10.0, 30.0
];

% Cholesky factor (lower triangular)
chol_sig = chol([1.0 0.5; 0.5 1.0], 'lower');

% Structure to hold data
WindDir.Data = wind_data;
WindDir.CholSig = chol_sig;

% Turbine indices and time
iT = [1, 2]';
t = 7.0;

phi = getWindDirT(WindDir, iT, t);

% Test
assert(isequal(size(phi), [2,1]));
assert(abs(phi(1) - 25.84752589085377) < 1e-6);
assert(abs(phi(2) - 25.140544918823198128) < 1e-6);

rmpath('./WindField/Direction_Interpolation_wErrorCov');

disp('All tests passed.');

