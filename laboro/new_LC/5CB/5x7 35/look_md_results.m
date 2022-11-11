clear;
close all;
M = [];
for i = 1 : 1
    [FILENAME, PATHNAME] = uigetfile(['.txt']);
    cd(PATHNAME);
    M = [M; dlmread([PATHNAME, FILENAME], ',')];
end
% file header
% t,temperature,ConvEKinTemp(ekin),sum_mom_xy,sum_vel[0],sum_vel[1],sum_vel[2],sum_vel_up[0],sum_vel_up[1],sum_vel_up[2],sum_vel_dw[0],sum_vel_dw[1],sum_vel_dw[2],sum_fup[0],sum_fup[1],sum_fup[2],sum_fdw[0],sum_fdw[1],sum_fdw[2]
cVx = 5;
cVy = 6;
cVz = 7;

cVxup = 8;
cVxdw = 11;
cFxup = 14;
cFxdw = 17;
natoms = 1330;
tstep = 2 * 1e-15; %m/s
mn = mean(M,1)


start1 = 1;
mean_Vx = 1000*cumsum(M(start1:end,cVx)/natoms)./((1:size(M,1)-start1+1)');
int_Vx_dt = 1000*tstep*1e9*cumsum((M(start1:end,cVx)-mn(cVx))/natoms);
int_Vy_dt = 1000*tstep*1e9*cumsum((M(start1:end,cVy)-mn(cVy))/natoms);
figure(12);
plot (int_Vx_dt, int_Vy_dt);




minvel_x = min(M(start1:end,cVx));
minvel_y = min(M(start1:end,cVy));

maxvel_x = max(M(start1:end,cVx));
maxvel_y = max(M(start1:end,cVy));

divx = 100;
divy = 100;

W = zeros(divx+1,divy+1);
for i = start1:size(M,1)
    sx = 1+fix(((M(i,cVx) - minvel_x) / (maxvel_x - minvel_x)) * divx);
	sy = 1+fix(((M(i,cVy) - minvel_y) / (maxvel_y - minvel_y)) * divy);
    W(sx,sy) = W(sx,sy) + 1.0;			
end

figure (15)
surf(W)


start2 = 2000;
figure(1);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), mean_Vx(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x)');xlabel('t, ps');ylabel('m/s');
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), int_Vx_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x) dt');xlabel('t, ps');ylabel('nm');


figure(2);
subplot(2,1,1);
plot (0.001*M(:,1), M(:,cFxup));title('sum F_x up'); xlabel('t, ps');
subplot(2,1,2);
plot (0.001*M(:,1), M(:,cFxdw));title('sum F_x down'); xlabel('t, ps');
 
figure(3);
subplot(2,1,1);
plot (0.001*M(:,1), cumsum(M(:,cFxup)));title('cumsum F_x up'); xlabel('t, ps');
subplot(2,1,2);
plot (0.001*M(:,1), cumsum(M(:,cFxdw)));title('cumsum F_x down'); xlabel('t, ps');
 

mean_fup = cumsum(M(:,cFxup))./((1:size(M,1))');
mean_fdown = cumsum(M(:,cFxdw))./((1:size(M,1))');
start2 = 50000;

figure(4);
subplot(2,1,1);
plot (0.001*M(1+start2:end,1), mean_fup(1+start2:end), [0.001*M(1+start2,1) 0.001*M(end,1)],[0 0]);title('mean F_x up'); xlabel('t, ps');ylabel('10^1^2 N/mol');
subplot(2,1,2);
plot (0.001*M(1+start2:end,1), mean_fdown(1+start2:end), [0.001*M(1+start2,1) 0.001*M(end,1)],[0 0]);title('mean F_x down'); xlabel('t, ps');ylabel('10^1^2 N/mol');
 










start1 = 1;
mean_t = cumsum(M(start1:end,cFxup)-M(start1:end,cFxdw))./((1:size(M,1)-start1+1)');

figure(5);
subplot(2,1,1);
start2 = 1000;
plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
 
subplot(2,1,2);
start2 = 50000;
plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 












start1 = 1;
mean_dFx = cumsum(M(start1:end,cFxup)-M(start1:end,cFxdw))./((1:size(M,1)-start1+1)');

figure(6);
subplot(3,1,1);
start2 = 10000;
plot (0.001*M(start1+start2:end,1), mean_dFx(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,2);
start2 = 100000;
plot (0.001*M(start1+start2:end,1), mean_dFx(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,3);
start2 = 200000;
plot (0.001*M(start1+start2:end,1), mean_dFx(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 


start1 = 1;
mean_dVx = 1000*cumsum(2*(M(start1:end,cVxup)-M(start1:end,cVxdw))/natoms)./((1:size(M,1)-start1+1)');
start1 = 1;
int_dVx_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVxup)-M(start1:end,cVxdw))/natoms);
figure(7);
subplot(2,1,1);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), mean_dVx(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x up - V_x down)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), int_dVx_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x up - V_x down) dt');xlabel('t, ps');ylabel('nm');
 

start1 = 1;
dFx = M(start1:end,cFxup)-M(start1:end,cFxdw);

len = size(dFx,1);
means = [];
N = 20
wind = fix(len / N);
for i = 1 : N
    m = mean(dFx(1 + (i-1)*wind:i*wind));
    means = [means m];
end
mat_ozh = mean(means)
D = sum((means - mat_ozh).^2)/(size(means, 2)-1);
sigma = sqrt(D / size(means, 2))

