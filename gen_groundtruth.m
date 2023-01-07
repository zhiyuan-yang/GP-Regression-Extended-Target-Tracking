function [groundtruth] = gen_groundtruth(model)
% This function generates ground truth for an extended target
% Input:         model
% Output:        groundtruth    A 2*model.simuTime matrix
% Author: Zhiyuan Yang
v = [10;0];
%v = [3;0];
ori = v/norm(v);
psi = atan2(ori(2),ori(1));
groundtruth(:,1) = model.initx(1:2,1);
for i=2:1:model.simuTime
    groundtruth(1:2,i) = groundtruth(1:2,i-1) + v*model.T;
    groundtruth(3,i) = psi;
end
end