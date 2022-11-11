clear;
close all;

[FILENAME, PATHNAME] = uigetfile(['.txt']);
cd(PATHNAME);

M = dlmread([PATHNAME, FILENAME], ',');
% file header
% t,temperature,ConvEKinTemp(ekin),sum_mom_xy,sum_vel[0],sum_vel[1],sum_vel[2],sum_vel_up[0],sum_vel_up[1],sum_vel_up[2],sum_vel_dw[0],sum_vel_dw[1],sum_vel_dw[2],sum_fup[0],sum_fup[1],sum_fup[2],sum_fdw[0],sum_fdw[1],sum_fdw[2]
cVx = 5;
cVx_wk_up = 8;
cVx_wk_dw = 11;
cVx_up = 14;
cVx_dw = 17;
cFx_wk_up = 20;
cFx_wk_dw = 23;
natoms_in_free_mols = 950;
natoms_worked = 10;
tstep = 2 * 1e-15; %m/s
mean(M,1)


% start1 = 1;
% mean_Vx = 1000*cumsum(M(start1:end,cVx)/natoms)./((1:size(M,1)-start1+1)');
% int_Vx_dt = 1000*tstep*1e9*cumsum(M(start1:end,cVx)/natoms);
% 
% start2 = 2000;
% figure(1);
% subplot(2,1,1);
% plot (0.001*M(start1+start2:end,1), mean_Vx(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x)');xlabel('t, ps');ylabel('m/s');
% subplot(2,1,2);
% plot (0.001*M(start1+start2:end,1), int_Vx_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x) dt');xlabel('t, ps');ylabel('nm');
% 

figure(2);
subplot(2,1,1);
plot (0.001*M(:,1), M(:,cFx_wk_up), [0.001*M(1,1) 0.001*M(end,1)],[0 0],'black');title('sum F_x _w_k up'); xlabel('t, ps');
subplot(2,1,2);
plot (0.001*M(:,1), M(:,cFx_wk_dw), [0.001*M(1,1) 0.001*M(end,1)],[0 0],'black');title('sum F_x _w_k down'); xlabel('t, ps');
 
figure(3);
subplot(2,1,1);
plot (0.001*M(:,1), cumsum(M(:,cFx_wk_up)), [0.001*M(1,1) 0.001*M(end,1)],[0 0],'black');title('cumsum F_x _w_k up'); xlabel('t, ps');
subplot(2,1,2);
plot (0.001*M(:,1), cumsum(M(:,cFx_wk_dw)), [0.001*M(1,1) 0.001*M(end,1)],[0 0],'black');title('cumsum F_x _w_k down'); xlabel('t, ps');
 

mean_f_wk_up = cumsum(M(:,cFx_wk_up))./((1:size(M,1))');
mean_f_wk_down = cumsum(M(:,cFx_wk_dw))./((1:size(M,1))');
start2 = 50000;

figure(4);
subplot(2,1,1);
plot (0.001*M(1+start2:end,1), mean_f_wk_up(1+start2:end), [0.001*M(1+start2,1) 0.001*M(end,1)],[0 0],'black');title('mean F_x _w_k up'); xlabel('t, ps');ylabel('10^1^2 N/mol');
subplot(2,1,2);
plot (0.001*M(1+start2:end,1), mean_f_wk_down(1+start2:end), [0.001*M(1+start2,1) 0.001*M(end,1)],[0 0],'black');title('mean F_x _w_k down'); xlabel('t, ps');ylabel('10^1^2 N/mol');
 










% start1 = 1;
% mean_t = cumsum(M(start1:end,cFxup)-M(start1:end,cFxdw))./((1:size(M,1)-start1+1)');
% 
% figure(5);
% subplot(2,1,1);
% start2 = 1000;
% plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
%  
%  
% subplot(2,1,2);
% start2 = 50000;
% plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 












start1 = 1;
mean_dF_wk_x = cumsum(M(start1:end,cFx_wk_up)-M(start1:end,cFx_wk_dw))./((1:size(M,1)-start1+1)');

figure(6);
subplot(3,1,1);
start2 = 1000;
plot (0.001*M(start1+start2:end,1), mean_dF_wk_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0],'black');title('mean (F_x _w_k up - F_x _w_k down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,2);
start2 = 10000;
plot (0.001*M(start1+start2:end,1), mean_dF_wk_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0],'black');title('mean (F_x _w_k up - F_x _w_k down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,3);
start2 = 100000;
plot (0.001*M(start1+start2:end,1), mean_dF_wk_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0],'black');title('mean (F_x _w_k up - F_x _w_k down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 




figure(7);
subplot(2,1,1);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), 1000*M(start1+start2:end,cVx_wk_up), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('V_x _w_k up');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), 1000*M(start1+start2:end,cVx_wk_dw), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('V_x _w_k down');xlabel('t, ps');ylabel('m/s');
 

start1 = 1;
mean_V_wk_x_up = 1000*cumsum(2*(M(start1:end,cVx_wk_up))/natoms_worked)./((1:size(M,1)-start1+1)');
start1 = 1;
int_V_wk_x_up_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_wk_up))/natoms_worked);



figure(8);
subplot(2,1,1);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), mean_V_wk_x_up(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x _w_k up)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), mean_V_wk_x_up(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x _w_k up) dt');xlabel('t, ps');ylabel('nm');
 

start1 = 1;
mean_V_wk_dw_x = 1000*cumsum(2*(M(start1:end,cVx_wk_dw))/natoms_worked)./((1:size(M,1)-start1+1)');
start1 = 1;
int_V_wk_x_dw_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_wk_dw))/natoms_worked);



figure(9);
subplot(2,1,1);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), mean_V_wk_dw_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x _w_k down)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), int_V_wk_x_dw_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x _w_k down) dt');xlabel('t, ps');ylabel('nm');
 


start1 = 1;
mean_dV_wk_x = 1000*cumsum(2*(M(start1:end,cVx_wk_up)-M(start1:end,cVx_wk_dw))/natoms_worked)./((1:size(M,1)-start1+1)');
start1 = 1;
int_dV_wk_x_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_wk_up)-M(start1:end,cVx_wk_dw))/natoms_worked);



figure(10);
subplot(2,1,1);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), mean_dV_wk_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x _w_k up - V_x _w_k down)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
start2 = 2000;
plot (0.001*M(start1+start2:end,1), int_dV_wk_x_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x _w_k up - V_x _w_k down) dt');xlabel('t, ps');ylabel('nm');
 

start1 = 1;
dFx = M(start1:end,cFx_wk_up)-M(start1:end,cFx_wk_dw);

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

