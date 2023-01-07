function [model] = gen_model()
% This function generates the model for extended target tracking
% Author: Zhiyuan Yang

model.simuTime = 50;    %%simulation timestep
model.T = 1;   %%sample time
model.motionmodel = 1;  
model.numBasisAngles = 10;
model.xf = linspace(0,2*pi,model.numBasisAngles); %%basis points
model.plot = true;  %%whether plot 
model.dimx = 2;     %%xc yc 
model.dimz = 2;
model.lambda = 20;  %%poisson measurement rate
model.shape = 2;                        %%shape1=rectangular shape2=circle
model.length = 2;
model.width = 1;    %%size of rectangular extent target
model.r = 1;        %%size of circle extent target

%% Gaussian Process model

model.mean = 0;
model.sigmar = 0.8;
model.sigmaf = 2;
model.l = pi/8;  %%length scale factor
basisAngleArray = transpose(linspace(0, 2*pi, model.numBasisAngles+1));
basisAngleArray(end) = [];
model.basisAngleArray = basisAngleArray;

%calculate K(uf,uf) since it will be used frequently.
cov_basis = compute_GP_covariance(basisAngleArray,basisAngleArray, ...
    model.sigmaf, model.sigmar,model.l);      
% using cholesky decomposition to calculate the inverse becasue it
% is more stable, more see[1].
cov_basis = cov_basis + 1e-6*eye(model.numBasisAngles);%to avoid numerical error
chol_cov_basis = chol(cov_basis);
inv_chol_cov_basis = inv(chol_cov_basis);
inv_cov_basis = inv_chol_cov_basis * inv_chol_cov_basis';
model.inv_cov_basis = inv_cov_basis;


%% motion model
% x_k+1 = F * x_k + W 
% x_k = [x_ck , f_k^T]^T
% W ~ N(0,Q)

%initialize
%initialize state
initxc = [1;1];
initpsi = 0;
initv = [0;0];
initpsid = 0;
initxf = model.mean*ones(model.numBasisAngles,1);
model.initx = [initxc;initpsi;initv;initpsid;initxf];
%initialize covariance matrix
initPxc = 10*eye(2);
initPv = 1*eye(2);
initPpsi = 1e-5;
initPpsid = 1e-5;
initPf = compute_GP_covariance(basisAngleArray,basisAngleArray,...
                               model.sigmaf,model.sigmar,model.l);
model.initP = blkdiag(initPxc,initPpsi,initPv,initPpsid,initPf);



F_k = kron([1,model.T;0,1],eye(3));
alpha = 0.001;    %%dynamic change of target extent set see 
F_f = exp(-alpha*model.T)*eye(model.numBasisAngles);
model.F = blkdiag(F_k,F_f);

sigma_q = 0.01;
sigma_q_psi = 0.0001;
Q_k = kron([model.T^3/3,model.T^2/2;model.T^2/2,model.T],[sigma_q^2,0,0;0,sigma_q^2,0;0,0,sigma_q_psi^2]);
Q_f = (1-exp(-2*alpha*model.T))* initPf;
model.Q = blkdiag(Q_k,Q_f);   %%motion noise--how much we trust the motion model

%% measurement model
model.R = 0.01*eye(2); %%measurement noise


end
