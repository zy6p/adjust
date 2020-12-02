%% task 1
% 问题1：设计矩阵与约束矩阵的关系，并通过程序验证

% read file
fileID = fullfile(cd,'CBW.txt');% cd is path
data = importdata(fileID);% import data
enum = 28;% n = 28
side = data(1:enum,1:3);% side numbers
L = data(1:enum,4);% L is observation
weight = data(1:enum,5);% weight
xpos = data(enum+1:end,2);% para x
ypos = data(enum+1:end,3);% para y

% Initialization
A = zeros(enum,16);
para = [xpos;ypos];% all para
sigma0 = 3;
P = diag(weight./sigma0).^2;

% RLESS
delta_X = para(side(:,2))-para(side(:,3));% \delta X_jk
delta_Y = para(side(:,2)+8)-para(side(:,3)+8);% \delta Y_jk
S_jk = sqrt(delta_X.^2+delta_Y.^2);% side length
l = L-S_jk;
A(sub2ind(size(A),side(:,1),side(:,2))) = -delta_X./S_jk;% X_j
A(sub2ind(size(A),side(:,1),side(:,2)+8)) = -delta_Y./S_jk;% Y_j
A(sub2ind(size(A),side(:,1),side(:,3))) = delta_X./S_jk;% X_k
A(sub2ind(size(A),side(:,1),side(:,3)+8)) = delta_Y./S_jk;% Y_k
N = A'*P*A;
c = A'*P*l;

% Minimization and internal constraints
K = [ones(1,8),zeros(1,8);...
    zeros(1,8),ones(1,8);...
    para(9:16)',para(1:8)'];
kesai = (N+K'*K)\c+(N+K'*K)\K'/(K/(N+K'*K)*K')*K/(N+K'*K)*c;
G = inv(N+K'*K)+(N+K'*K)\K'/(K/(N+K'*K)*K')*K/(N+K'*K);

eee = A*kesai-L;

% conpare A and K
rkK = rank(K)
rkA = rank(N)
rank([N,K'])
fprintf("rk(K) = %d = m-q\nrk[A',K'] = rk(A)+rk(K) = %d+%d =  m\n设计矩阵与约束矩阵组成的矩阵满秩\n",rkK,rkA,rkK)
%% task 2
% 问题2：题中给的约束是否是内约束，如果不是如何由它得到内约束

K2 = (eye(16)-pinv(N)*N)*rand(16,3);
K2 = K2';
A*K'
fprintf("K不是内约束，A*K'~=0")
rank(K2)
A*K2'
if max(abs(A*K2'))<1e-10
    fprintf('K2是内约束\n求内约束方法：')
end
%% 
% $K^{\prime } =\left(I_m -N^- N\right)\alpha$  α为任意16*3的矩阵
%% task 3
% 问题3:用内约束解算最小二乘结果，对比minoless结果，以及对比两者的协因数阵

kesai2 = (N+K2'*K2)\c+(N+K2'*K2)\K2'/(K2/(N+K2'*K2)*K2')*K2/(N+K2'*K2)*c;
G2 = inv(N+K2'*K2)+(N+K2'*K2)\K2'/(K2/(N+K2'*K2)*K2')*K2/(N+K2'*K2);
Q2 = (N+K2'*K2)\N/(N+K2'*K2);
kesai_minoless = N*pinv(N*N)*c;
Q_minoless = pinv(N);
if Q_minoless-Q2<1e-10
    fprintf('两者的协因数阵相等')
end
%% task 4
% 问题4:通过S转换，使任意一个解（N的任意广义逆与c的乘积）得到minoless的解

[U,lambda,V] = eig(N);
delta = diag(lambda(abs(lambda)>1e-10))^-1;
N_inv = U*[delta,rand(13,3);rand(3,13),rand(3)]*V';
kesai3 = N_inv*c;
kesai4 = pinv(N)*N*kesai3;
if kesai4-kesai_minoless<1e-10
    fprintf('通过S转换成功')
end
%% task 5
% 问题5：Minoless的解为 N（NN）-c，计算是否和pinv（N）*c相等。如果是，那么N（NN）-等于pinv（N）么，为什么？

kesai5 = pinv(N)*c;
if kesai5-kesai_minoless<1e-10
    fprintf('两者的解相等')
end
[U2,lambda2,V2]  =eig(N*N);
delta2 = diag(lambda2(abs(lambda2)>1e-10))^-1;
N_minoless = U2*[delta2,rand(13,3);rand(3,13),rand(3)]*V2';
if abs(N_minoless-pinv(N))>1e-10
    fprintf('N（NN）- 和pinv（N）不相等')
end
%% 
% 因为N（NN）-和pinv（N）都是$N_{\mathrm{rs}}^-$，而对称反射广义逆不唯一
%% task 6
% 问题6 （选做）：以上都是用最小二乘求解，请使用Huber函数求解，将两个结果的残差放到一个表格对比

c = 3;
beta = kesai2;
W = P;

e = A*beta-L;
sigma = (e'*W*e./(size(L,1)-2))^0.5;
x_e = e/sigma;% x in Huber-M
W = Huber_M_Weight(x_e,P, c);% W in Huber-M
% W can also be changed one by one
delta_beta = -(A'*W*A)\A'*W*e;
beta = beta+delta_beta;

ee = [eee,e]'
emean = mean(ee);
ee = ee-emean;
plot(side(:,1),ee(1,:)','d:')
hold on
plot(side(:,1),ee(2,:)','^--')
hold off

legend('show')
%%
function weight = Huber_M_Weight(x_e,P, c)
x_e = abs(x_e);
weight = zeros(size(P));
for ii = 1:size(P,1)
    weight(ii,ii) = P(ii,ii)*(x_e(ii)<c)+P(ii,ii)/x_e(ii)*(x_e(ii)>=c);
end
end