%     dir_mode = Direction_Constant()
%     WindDir = 270
%     iT = [1, 2, 3]
% %     phi = getWindDirT(dir_mode, WindDir, iT, nothing)
%     for ph in phi
%         @test ph ≈ 270.0
%     end

addpath('./WindField/Direction_Constant');
WindDir = 270;
iT = [1, 2, 3];
phi = getWindDirT(WindDir, iT);

tol = 1e-8; % set a tolerance for floating-point comparison
for i = 1:length(phi)
    ph = phi(i);
    assert(abs(ph - 270.0) < tol);
end


