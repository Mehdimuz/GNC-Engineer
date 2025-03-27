clear; close all; clc
MC.run = true;

%%% Addmissible Variations
% Add whatever you want to this to  randomize every run
% 
MC.q = [-1 1];
MC.w = [-.01 .01];
MC.J = false;

MC.runs = 1;
while MC.runs < 40
    mission_simulation

    % Make a log that persists through every run
    LOGLOG.q(:,:,MC.runs) = LOG.trueq;
    LOGLOG.q_err(:,:,MC.runs) = LOG.qerr;
    LOGLOG.q_err_est(:,:,MC.runs) = LOG.q_err_est;
    LOGLOG.q_err_con(:,:,MC.runs) = LOG.q_err_con;

    LOGLOG.sig3(:,:,MC.runs) = LOG.sig3;

    LOGLOG.h(:,:,MC.runs) = LOG.h;

    MC.runs = MC.runs + 1;
end
%% 
%%%%% plot some 
% just a few examples
% these are 3-D arrays now so it can be a little difficult

% overall attitude estimation error over time
figure(1)
hold off, plot(tspan/60,2*acosd( permute(abs(LOGLOG.q_err(4,:,:)) , [2 3 1]) ),'HandleVisibility','off')
hold on, plot(tspan/60, 180/pi*vecnorm( mean(LOGLOG.sig3(1:3,:,:),3) ), '--k')
legend("Average of all 3 \sigma bounds"), grid on
title("Monte Carlo Attitude Estimate Error"), xlabel("Time [m]"), ylabel("Error [°]")

% total attitude control error
% figure(1)
% hold off, plot(tspan/60,2*acosd( permute(abs(LOGLOG.q_err_con(4,:,:)) , [2 3 1]) ))
% hold on, plot(xlim, [6.3 6.3], '--k')

% wheel saturation
figure(2)
subplot(311)
plot(tspan/60, permute(LOGLOG.h(1,:,:) , [2 3 1])/POLAR.wheel_J(1)/POLAR.maxw_wheel(1) )
ylim([-1 1]), grid on
title("Wheel Saturation Over Time")
subplot(312)
plot(tspan/60, permute(LOGLOG.h(2,:,:) , [2 3 1])/POLAR.wheel_J(2)/POLAR.maxw_wheel(2) )
ylim([-1 1]), grid on
subplot(313)
plot(tspan/60, permute(LOGLOG.h(3,:,:) , [2 3 1])/POLAR.wheel_J(3)/POLAR.maxw_wheel(3) )
ylim([-1 1]), grid on
xlabel("Time [minutes]")

