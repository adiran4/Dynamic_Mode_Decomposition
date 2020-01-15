clear all; clc
data_matrix = [];
v = VideoReader('IMG_0831.MOV');
vidFrames = read(v);
numFrames = get(v, 'numberOfFrames');
for i = 1:numFrames
frame = vidFrames(:, :, :, i);
frame = double(rgb2gray(frame));
data_matrix(i,:) = reshape(imresize(frame, [264, 400]),[1,264*400]);
end
%%
data_matrix = data_matrix.';
X1 = data_matrix(:, 1:end - 1);
X2 = data_matrix(:, 2:end);
[U, S, V] = svd(X1, 'econ');
%%
semilogy(diag(S))
xlabel('Singular values')
ylabel('Amplitude(log scale)')
title('Singular Value Spectrum')
%%
r = 4;
U = U(:, 1:r);
S = S(1:r, 1:r);
V = V(:, 1:r);
Atlide = U'*X2*V/S;
[W, D] = eig(Atlide);
Phi = X2*V/S*W;
%%
dt = 1/3.9;
mu = diag(D);
omega = log(mu)/dt;
%%
x1 = data_matrix(:, 1);
b = Phi\x1;
mm1 = size(X1, 2);
u_modes = zeros(r, mm1);
t = (0:mm1 - 1)*dt;
for iter = 1:mm1
u_modes(:, iter) = (b.*exp(omega*t(iter)));
end
X_dmd = Phi*u_modes;
%%
for k = 1:225
y = reshape(X_dmd(:, k), [264 400]);
imshow(uint8(y));
pause(0.001);
end
%%
X_sparse = X1 - abs(X_dmd);
for m = 1:225
y = reshape(X_sparse(:, m), [264, 400]);
imshow(uint8(y));
pause(0.001);
end
%%
R_Matrix = X_sparse.*(X_sparse < 0);
X_sparse_dmd = X_sparse - R_Matrix;
X_bg = R_Matrix + abs(X_dmd);
X_fg = X_sparse - R_Matrix;
for m = 1:225
y = reshape(X_sparse_dmd(:, m), [264, 400]);
imshow(uint8(y));
pause(0.001);
end
