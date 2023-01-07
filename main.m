%  Code for Gaussian Process Regression Extended Target Tracking
%  Author: Zhiyuan Yang
%  Date: 2022.12.02

%%
%close all;
%clear all;
%clc;

model = gen_model();                    %%generate model
truth = gen_groundtruth(model);   %%generate groundtruth
meas = gen_meas(model,truth);     %%generate measurements
model.plot = false;

%% run the filter
est = gp_ekf_filter(meas,model);

numBasisAngles = model.numBasisAngles;
basisAngleArray = model.basisAngleArray;
model.plot = true;
l = model.length;
w = model.width;
if model.plot
    figure(1);
    hold on;
    for k = 1:1:model.simuTime
        p1 = plot(truth(1,k),truth(2,k),'ro');
        p2 = plot(est.x(1,k),est.x(2,k),'gs');
        psi = truth(3,k);
        M = [cos(psi),-sin(psi);sin(psi),cos(psi)];
        if model.shape == 1
            %plot rectangular extent
            edg1 = [linspace(-l/2,l/2,100);w/2*ones(1,100)];
            edg1 = M*edg1 + truth(1:2,k);
            p3 = plot(edg1(1,:),edg1(2,:),'-k','LineWidth',2);

            edg2 = [linspace(-l/2,l/2,100);-w/2*ones(1,100)];
            edg2 = M*edg2 + truth(1:2,k);
            plot(edg2(1,:),edg2(2,:),'-k','LineWidth',2);

            edg3 = [-l/2*ones(1,100);linspace(-w/2,w/2,100)];
            edg3 = M*edg3 + truth(1:2,k);
            plot(edg3(1,:),edg3(2,:),'-k','LineWidth',2);

            edg4 = [l/2*ones(1,100);linspace(-w/2,w/2,100)];
            edg4 = M*edg4 + truth(1:2,k);
            plot(edg4(1,:),edg4(2,:),'-k','LineWidth',2);

               
        elseif model.shape ==2
            angle = linspace(0,2*pi,300);
            p3 = plot(truth(1,k)+model.r*cos(angle),truth(2,k)+model.r*sin(angle),'-k','LineWidth',2);
        end
        %plot measurements
        for i=1:1:length(meas{k})
            p4 = plot(meas{k}(1,i),meas{k}(2,i),'bx');
        end
        %plot predicted extent
        xxf = est.x(7:end,k) .* cos(basisAngleArray-est.x(3,k)) + est.x(1,k);
        yxf = est.x(7:end,k) .* sin(basisAngleArray-est.x(3,k)) + est.x(2,k); 
        p5 = plot(xxf,yxf,'m-','Linewidth',2);
    end
    legend([p1,p2,p3,p4,p5],...
        {'true position','predicted position','true extent','measurements','predicted extent'});
    xlabel('X/m');
    ylabel('Y/m');
    axis equal;
    hold off;
end

