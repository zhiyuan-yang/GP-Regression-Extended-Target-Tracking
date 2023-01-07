function [meas] = gen_meas(model, groundtruth)
%Input: groundtruth     A dimx*model.simuTime matrix
%       model
%Output: meas           A simuTime*1 cell. Each cell meas{k} contains a
%                       dimz*num_meas matrix of measurements at time k. 
%Author: Zhiyuan Yang


rng("shuffle");  %%rand number
w = model.width;
l = model.length;
theta_rec = atan2(w,l);
meas = cell(model.simuTime, 1);
for k=1:1:model.simuTime
    num_meas_k = poissrnd(model.lambda); %%number of measurements
    theta_k = 2*pi*rand(1,num_meas_k);
    psi = groundtruth(3,k);
    M = [cos(psi),-sin(psi);sin(psi),cos(psi)];  %%rotation matrix
    if model.shape == 1
        theta_k = theta_k - pi;
        id1 = find(theta_k>theta_rec & theta_k<pi-theta_rec); 
        z1 =  groundtruth(1:2,k) + w/2*[1./tan(theta_k(id1));ones(1,length(id1))] + 0.01*rand(2,length(id1));
        id2 = find(theta_k<-theta_rec &theta_k>-pi+theta_rec);
        z2 =  groundtruth(1:2,k) + w/2*[-1./tan(theta_k(id2));-ones(1,length(id2))] + 0.01*rand(2,length(id2));
        id3 = find(theta_k<theta_rec & theta_k>-theta_rec);
        z3 = groundtruth(1:2,k) + l/2*[ones(1,length(id3));tan(theta_k(id3))]+0.01*rand(2,length(id3));
        id4 = find(theta_k<-pi+theta_rec | theta_k>pi-theta_rec);
        z4 = groundtruth(1:2,k) + l/2*[-ones(1,length(id4));tan(theta_k(id4))]+0.01*rand(2,length(id4));
        meas{k} = [z1,z2,z3,z4];
    elseif model.shape == 2
        z = groundtruth(1:2,k) + model.r*[cos(theta_k) ; sin(theta_k)] + 0.1*rand(2,num_meas_k);  
        meas{k} = [meas{k},z];
    end
    meas{k} = M * (meas{k}-groundtruth(1:2,k)) + groundtruth(1:2,k);
end
