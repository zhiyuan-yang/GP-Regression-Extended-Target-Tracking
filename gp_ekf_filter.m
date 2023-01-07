function [est] = gp_ekf_filter(meas,model)
% This function runs EKF filter with GP regression in k scan.
% Input:    meas       A model.simuTime*1 cell
%           model      Struct
% Output:   est        A struct
%           est.predx  A dimx_k*simuTime matrix
%           est.predP  A dimx_k*dimx_k*simuTime matrix
% stateVector    x_k = [x_ck; y_ck; psi_k; x_ck'; y_ck'; psi_k'; f_k]
%                f_k = [f0_k; f1_k;...; fn_k]
% motion model                x_k+1 = F*x_k + w_k      w_k~N(0,Q)
% measurement model           z_k = h_k(x_k) + e_k     e_k~N(0,R)
% Reference: [1] N. Wahlström and E. Özkan, ‘Extended Target Tracking Using Gaussian Processes’, 
%                IEEE Transactions on Signal Processing, vol. 63, no. 16,
%                pp. 4165–4178, Aug. 2015.
% Author: Zhiyuan Yang

F = model.F;
Q = model.Q;
R = model.R;
sigmaf = model.sigmaf;
sigmar = model.sigmar;
l = model.l;
dimx = model.dimx;
dimz = model.dimz;
numBasisAngles  = model.numBasisAngles;
basisAngleArray = model.basisAngleArray;
inv_cov_basis = model.inv_cov_basis;
est.prex(:,1) = model.initx;
est.x(:,1)= model.initx;
est.preP(:,:,1) = model.initP;
est.P(:,:,1) = model.initP;

for k = 2:1:model.simuTime
%Extended Kalman Filter
%%Prediction
    est.prex(:,k)= F * est.x(:,k-1);
    est.preP(:,:,k) = F * est.P(:,:,k-1) * F' + Q;
    
    xc = est.prex(1:dimx,k);
    xpsi = est.prex(dimx+1,k);
    xf = est.prex(2*(dimx+1)+1:end,k);

    Hk = zeros(dimz*length(meas{k}),length(est.prex(:,k)));
    Rk = zeros(dimz*length(meas{k}),dimz*length(meas{k}));
    prez = zeros(dimz,length(meas{k}));
%%Update
    for i=1:1:length(meas{k}) %go through all measurements
        %cal predictz
        
        diffvector = meas{k}(:,i)-xc;
        orientation_z = diffvector./norm(diffvector);
        globalangle = atan2(orientation_z(2),orientation_z(1));
        localangle = mod((globalangle - xpsi),2*pi);
        cov_meas_basis = compute_GP_covariance(localangle,basisAngleArray,sigmaf,sigmar,l);   
        
        Hf = cov_meas_basis * inv_cov_basis;
        prez(:,i) = xc + orientation_z * Hf * xf;
        
        %calculate Rk
        meas_cov = compute_GP_covariance(localangle,localangle,sigmaf,sigmar,l);
        R_i = meas_cov - cov_meas_basis*inv_cov_basis*cov_meas_basis';
        Rk(1+(i-1)*dimz:i*dimz,1+(i-1)*dimz:i*dimz) = R^2*eye(2)+orientation_z...
            *R_i*orientation_z';

        %linearize the model see Appendix B
        dHf_du= -1/l^2*sin(localangle-basisAngleArray').*cov_meas_basis*inv_cov_basis; %(57c-57e)
        dp_dw=  (diffvector*diffvector')./norm(diffvector)^3 - 1/norm(diffvector)*eye(dimx);%(57b)
        dtheta_dw= 1/norm(diffvector)^2*[diffvector(2),-diffvector(1)];%(57a)

        dh_dx_c= eye(2) + dp_dw*(Hf*xf) + orientation_z*dtheta_dw*(dHf_du*xf);%(56a)
        dh_dx_psi=-orientation_z*dHf_du*xf;%(56b)
        dh_dx_f =orientation_z*Hf;%(56c)

        Hk(1+(i-1)*dimz:i*dimz,:)= [dh_dx_c,dh_dx_psi,zeros(2,3),dh_dx_f]; 
        % Hk is a (dimz*nummeas)*dimxk matrix 
    end
    currmeas = meas{k}(:);       %form[x1;y1;x2;y2...xn;yn] 2l*1 vector
    prez = prez(:);
    
    Sk = Hk*est.preP(:,:,k)*Hk' + Rk;
    Kk = est.preP(:,:,k) * Hk'/Sk;
    est.x(:,k) = est.prex(:,k) + Kk * (currmeas - prez);
    est.x(3,k) = 0;
    est.P(:,:,k) = (eye(length(est.prex(:,k))) - Kk*Hk)*est.preP(:,:,k);
    est.P(:,:,k) = (est.P(:,:,k)+est.P(:,:,k)')/2;
end

end