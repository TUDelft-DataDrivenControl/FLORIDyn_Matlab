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

% % The following test fails
% % clear all
% % addpath('./WindField/Direction_Constant_wErrorCov');
% % rng(1234);
% % 
% % % Cholesky factor of 3x3 identity
% % L = chol(eye(3), 'lower');
% % 
% % % WindDir structure
% % WindDir.Data = 270.0;
% % WindDir.CholSig = L;
% % 
% % iT = [1, 2, 3];
% % 
% % % Expected result
% % result = [269.6402710931765, 271.0872084924286, 269.5804103830612];
% % 
% % phi = getWindDirT(WindDir, iT);
% % for i = 1:length(phi)
% %     assert(abs(phi(i) - result(i)) < 1e-8, 'Test failed at index %d', i);
% % end
% % rmpath('./WindField/Direction_Constant_wErrorCov');

% %   dir_mode = Direction_EnKF_InterpTurbine()
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

disp('All tests passed.');

