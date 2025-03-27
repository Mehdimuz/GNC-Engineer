function SAT = Wheel_Saturation(SAT, L, dt)

% SAT is the satellite struct
% L is the desired effective torque. A positive effective torque on the
% satellite is a negative torque on the wheels
MAXW = SAT.maxw_wheel; MINW = -SAT.maxw_wheel;
MAXT = SAT.maxT_wheel; MINT = -SAT.maxT_wheel;

wheeltorques = min(L, MAXT);
wheeltorques = max(wheeltorques, MINT);
wheeltorques = min(wheeltorques,(MAXW-SAT.wheelspeeds).*SAT.wheel_J/dt);
wheeltorques = max(wheeltorques,(SAT.wheelspeeds-MAXW).*SAT.wheel_J/dt);

SAT.wheeltorques = wheeltorques;
end

