close all;
c = {"Roll","Pitch","Yaw"};
for i = 1:3
subplot(3,1,i)
plot(tspan/60,LOG.qerr(4,:).*asin(LOG.qerr(i,:))*2*180/pi,'linewidth', 1)       % Estimation Error Plots
ylabel(c{i} +" (Degrees)")
hold on
plot(tspan/60,[-1; 1]*LOG.sig3(i,:)*180/pi, '--r', "linewidth", 2)              % 3-sigma bounds 
ylim(180/pi*[-3 3]*(LOG.sig3(i,end)))
xticks([0:floor(max(tspan)/60)])
grid on
legend("Attitude Error","3-\sigma bounds")
end
xlabel("Time (Minutes)")
subplot(3,1,1), title("Attitude Estimation Error",'FontSize',18)

figure(2)
c = {"Boresight","Spin Axis"};
for i = 1:2
subplot(2,1,i)
plot(tspan/60,(asind(LOG.q_err_con(i,:).*LOG.q_err_con(4,:))*2),'linewidth', 2)   % Pointing/Ang Rate Error Plots
ylabel(c{i} +" Tracking Error (Degrees)", "FontSize",14)
%ylim(max(abs(ylim))*[-1 1])
ylim([-6.3 6.3])
xticks([0:floor(max(tspan)/60)])
grid on
hold on
plot(tspan/60,(LOG.w_des(i,:)-LOG.omega(i,:))*180/pi,'linewidth', 0.8,'Color',[0.6 0.6 0.6])
plot(xlim', 6.2*[1;1]*[-1 1], '--r','linewidth',2)
legend("Tracking Error","Angular Rate Error","Pointing Requirement",'fontsize',13)
ylim([-10 10])
end
xlabel("Time (Minutes)")
subplot(2,1,1), title("Pointing error",'FontSize',18)


%% TARGETER CHECK
%while doAnimation == true

if doAnimation
figure, pause(.1)
set(gcf,"Position",[100 100 1000 700])
subplot(1,3,1:2)
[X,Y,Z] = sphere(60);
earth = surf(X*6371, Y*6371, Z*6371, "EdgeAlpha",0, 'FaceAlpha',.6,...
    'FaceColor','texturemap','CData',flip(imread('Earth Texture.jpg')));

hold on
RR = 1;
r = LOG.position(:,RR); v = LOG.velocity(:,RR);
los = LOG.los(:,RR);    los_dot = LOG.los_d(:,RR);
pol = surf(X*60+r(1),Y*60+r(2),Z*60+r(3),...
    'FaceColor',[.8,.8,.8], 'EdgeAlpha',0);
obj = surf(X*30+r(1)+los(1), Y*30+r(2)+los(2), Z*30+r(3)+los(3),...
    'FaceColor','red','EdgeAlpha',0);

dir = plot3(r(1)+los(1)+[0 los_dot(1)]*300, ...
      r(2)+los(2)+[0 los_dot(2)]*300, ...
      r(3)+los(3)+[0 los_dot(3)]*300,'--k','LineWidth',2);

A = q2a(LOG.q_des(:,RR));
fx = plot3(r(1)+A(1,1)*[0,400],r(2)+A(1,2)*[0, 400],r(3)+A(1,3)*[0,400], 'r','linewidth',2);
fy = plot3(r(1)+A(2,1)*[0,400],r(2)+A(2,2)*[0, 400],r(3)+A(2,3)*[0,400], 'g','linewidth',2);
fz = plot3(r(1)+A(3,1)*[0,400],r(2)+A(3,2)*[0, 400],r(3)+A(3,3)*[0,400], 'b','linewidth',2);

lims = max(vecnorm(LOG.los));
xlim(r(1)+2*[-1 1]*lims)   % PDR EDIT 
ylim(r(2)+2*[-1 1]*lims)
zlim(r(3)+2*[-1 1]*lims)
% axis equal

hold off, axis off, set(gcf, 'color',[.5 .5 .5])



subplot(1,3,3)
object_trail = polarplot(nan(50,1), nan(50,1)); hold on
object_mark = polarplot(nan, nan,'x', 'MarkerSize',20);
stars_plot = polarplot(nan,nan,'.');
rlim([0 6.2]);
rticks([2 4 6]);
thetaticks([0 90 180 270]), hold off
set(gca, 'color',[.1 .1 .1])
set(gca, 'GridColor',[.9 .9 .9])
set(gca, 'RColor', [.7 .7 .7])

subplot(1,3,1:2)
for RR = 1:Lsim
    r = LOG.position(:,RR); v = LOG.velocity(:,RR);
    %los = LOG.los(:,RR);    los_dot = LOG.los_d(:,RR);  %ORIGINAL
    los = LOG.los(:,RR);    los_dot = LOG.los_d(:,RR); % PDR EDIT --> this way appears to eventually settle on the object... 

    pol.XData = X*60+r(1);
    pol.YData = Y*60+r(2);
    pol.ZData = Z*60+r(3);

    obj.XData = X*30+r(1)+los(1);
    obj.YData = Y*30+r(2)+los(2);
    obj.ZData = Z*30+r(3)+los(3);
    
    dir.XData = r(1)+los(1)+[0 los_dot(1)]*300;
    dir.YData = r(2)+los(2)+[0 los_dot(2)]*300;
    dir.ZData = r(3)+los(3)+[0 los_dot(3)]*300;
    
    A = q2a(LOG.trueq(:,RR));
    fx.XData = r(1)+A(1,1)*[0,400];
    fx.YData = r(2)+A(1,2)*[0,400];
    fx.ZData = r(3)+A(1,3)*[0,400];
    fy.XData = r(1)+A(2,1)*[0,400];
    fy.YData = r(2)+A(2,2)*[0,400];
    fy.ZData = r(3)+A(2,3)*[0,400];
    fz.XData = r(1)+A(3,1)*[0,400];
    fz.YData = r(2)+A(3,2)*[0,400];
    fz.ZData = r(3)+A(3,3)*[0,400];

    xlim(r(1)+2*[-1 1]*lims)   % PDR EDIT 
    ylim(r(2)+2*[-1 1]*lims)
    zlim(r(3)+2*[-1 1]*lims)
    view(atan2d(r(2),r(1))+10, 20)

    % xlim(r(1)+[-500 500])
    % ylim(r(2)+[-500 500])
    % zlim(r(3)+[-500 500])
    drawnow limitrate
    title(sprintf("%d:%02d",floor(RR*FC.dt/60),mod(floor(RR*FC.dt),60)))

    b(:,RR) = A*los;
    bd(:,RR) = A*los_dot;

    err_mag = acosd(dot([0;0;1], b(:,RR))./norm(b(:,RR)));
    err_the = atan2(b(1,RR),-b(2,RR));
    object_mark.RData = err_mag; object_mark.ThetaData = err_the;
    object_trail.RData = [object_trail.RData(2:end), err_mag];
    object_trail.ThetaData = [object_trail.ThetaData(2:end), err_the];

    bstars = A*stars;
    inds = find( ([0 0 1]*bstars) >= cosd(6.2));
    stars_plot.RData = acosd([0 0 1]*bstars(:,inds));
    stars_plot.ThetaData = atan2(bstars(1,inds), -bstars(2,inds));

end

end

%end
%% Power

figure('Position',[100 100 1000 600])
subplot(211)
plot(tspan/3600, [LOG.pwr_gen; LOG.pwr_con],'linewidth',2)
ylim([0 max(ylim*1.2)]), grid on
xlabel("Time (hours)"), ylabel("Wattage")
title(["Power Usage in " + FC.mode + " Mode"])
legend("Generated","Consumed","Location","best")
subplot(212)
plot(tspan/3600, LOG.bat, 'linewidth', 2)
xlabel("Time (hours)"), ylabel("Battery Charge (W hr)")
ylim([0 PWR.bat_charge_max]), grid on
%% Ground Track
figure
imshow(imread('Earth Texture.jpg'))
hold on
plot((LOG.lon(~LOG.TX)+180)*1022/360, -(LOG.lat(~LOG.TX)-90)*511/180,...
    '.r', 'MarkerSize', 15,'HandleVisibility','off')
hold on
plot((LOG.lon(LOG.TX)+180)*1022/360, -(LOG.lat(LOG.TX)-90)*511/180,'.g', 'MarkerSize', 15)
legend("Transmitting",'Fontsize',15)
title("Ground Track Showing Location When Transmitting", 'FontSize',20)
