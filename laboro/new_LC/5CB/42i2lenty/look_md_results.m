fm = menu('Begin?', 'Yes', 'No');
if (fm ~= 1)
    return;
end


clear;
close all;
M = [];
for i = 1 : 8
    [FILENAME, PATHNAME] = uigetfile(['.txt']);
    if ((0 == FILENAME) & (0 == PATHNAME))
        break;
    end
    cd(PATHNAME);FILENAME
    M = [M; dlmread([PATHNAME, FILENAME], ',')];
end
% file header
% t,temperature,ConvEKinTemp(ekin),sum_mom_xy,sum_vel[0],sum_vel[1],sum_vel[2],sum_vel_up[0],sum_vel_up[1],sum_vel_up[2],sum_vel_dw[0],sum_vel_dw[1],sum_vel_dw[2],sum_fup[0],sum_fup[1],sum_fup[2],sum_fdw[0],sum_fdw[1],sum_fdw[2]
cVx = 5;
cVy = 6;
cVz = 7;

ct = 1;

cVx_wk_up = 8;
cVx_wk_dw = 11;
cVx_up = 14;
cVx_dw = 17;
cFx_wk_up = 20;
cFx_wk_dw = 23;
natoms = 1932;
natoms_worked = 12;
tstep = 2 * 1e-15; %m/s
mean(M,1)

DT = 100;
start1 = 1000;
t = M(start1:DT:end,ct);
v_Vx = M(start1:DT:end,cVx);
v_Vx_wk_dw = M(start1:DT:end,cVx_wk_dw);
v_Vx_wk_up = M(start1:DT:end,cVx_wk_up);

m_Vx = mean(v_Vx);
m_Vx_wk_dw = mean(v_Vx_wk_dw);
m_Vx_wk_up = mean(v_Vx_wk_up);

mv_Vx = v_Vx - m_Vx;
mv_Vx_wk_dw = v_Vx_wk_dw - m_Vx_wk_dw;
mv_Vx_wk_up = v_Vx_wk_up - m_Vx_wk_up;

[K,tau]=correlation_function(t, mv_Vx.*mv_Vx, mv_Vx_wk_dw);
figure(27)
plot(tau, K);

[K,tau]=correlation_function(t, mv_Vx.*mv_Vx, mv_Vx_wk_up);
figure(28)
plot(tau, K);

[K,tau]=correlation_function(t, mv_Vx_wk_dw, mv_Vx_wk_up);
figure(29)
plot(tau, K);

% len = 10000*fix(size(mv_Vx, 1)/10000)
% 
% kf = zeros(len,1);
% kf_wk_dw = zeros(len,1);
% kf_wk_up = zeros(len,1);

% for i = 1 : len    
%     kf(i) = korrelation_momentum(i-1, mv_Vx, mv_Vx_wk_dw);    
% end




start1 = 1;
mean_Vx = 1000*cumsum(M(start1:end,cVx)/natoms)./((1:size(M,1)-start1+1)');

int_Vx_dt = 1000*tstep*1e9*cumsum(M(start1:end,cVx)/natoms);
int_Vy_dt = 1000*tstep*1e9*cumsum(M(start1:end,cVy)/natoms);
figure(2);
plot (int_Vx_dt, int_Vy_dt);


start2 = 1000;




min_Vx = min(M(start1+start2:end,cVx));
min_Vx_wk_dw = min(M(start1+start2:end,cVx_wk_dw));
min_Vx_wk_up = min(M(start1+start2:end,cVx_wk_up));

max_Vx = max(M(start1+start2:end,cVx));
max_Vx_wk_dw = max(M(start1+start2:end,cVx_wk_dw));
max_Vx_wk_up = max(M(start1+start2:end,cVx_wk_up));

div_Vx = 199;
div_Vx_wk_dw = 99;
div_Vx_wk_up = 99;

 
W_dw = zeros(div_Vx_wk_dw+1,div_Vx+1);
W_up = zeros(div_Vx_wk_up+1,div_Vx+1);
for i = start1+start2:size(M,1)
    s_Vx = 1+fix(((M(i,cVx) - min_Vx) / (max_Vx - min_Vx)) * div_Vx);
	s_wk_dw = 1+fix(((M(i,cVx_wk_dw) - min_Vx_wk_dw) / (max_Vx_wk_dw - min_Vx_wk_dw)) * div_Vx_wk_dw);
	s_wk_up = 1+fix(((M(i,cVx_wk_up) - min_Vx_wk_up) / (max_Vx_wk_up - min_Vx_wk_up)) * div_Vx_wk_up);
    W_dw(s_wk_dw,s_Vx) = W_dw(s_wk_dw,s_Vx) + 1.0;			
    W_up(s_wk_up,s_Vx) = W_up(s_wk_up,s_Vx) + 1.0;			
end

s_Vx_0 = 1+fix(((0.0 - min_Vx) / (max_Vx - min_Vx)) * div_Vx);
s_wk_dw_0 = 1+fix(((0.0 - min_Vx_wk_dw) / (max_Vx_wk_dw - min_Vx_wk_dw)) * div_Vx_wk_dw);
s_wk_up_0 = 1+fix(((0.0 - min_Vx_wk_up) / (max_Vx_wk_up - min_Vx_wk_up)) * div_Vx_wk_up);

% figure (23)
% surf(W)

figure (24)
contour(W_dw)
hold on;
plot ([s_Vx_0 s_Vx_0], [0 div_Vx_wk_dw+1], '-');
plot ([0 div_Vx+1 ], [s_wk_dw_0 s_wk_dw_0], '-');
hold off;

figure (25)
contour(W_up)
hold on;
plot ([s_Vx_0 s_Vx_0], [0 div_Vx_wk_up+1], '-');
plot ([0 div_Vx+1 ], [s_wk_up_0 s_wk_up_0], '-');
hold off;

% figure(21);
% plot (M(start1+start2:end,cVx), M(start1+start2:end,cVx_wk_dw), '.');
% 
% figure(22);
% plot (M(start1+start2:end,cVx), M(start1+start2:end,cVx_wk_up), '.');




% 
% minvel_x = min(M(start1:end,cVx));
% minvel_y = min(M(start1:end,cVy));
% 
% maxvel_x = max(M(start1:end,cVx));
% maxvel_y = max(M(start1:end,cVy));
% 
% divx = 100;
% divy = 100;
% 
% W = zeros(divx+1,divy+1);
% for i = start1:size(M,1)
%     sx = 1+fix(((M(i,cVx) - minvel_x) / (maxvel_x - minvel_x)) * divx);
% 	sy = 1+fix(((M(i,cVy) - minvel_y) / (maxvel_y - minvel_y)) * divy);
%     W(sx,sy) = W(sx,sy) + 1.0;			
% end
% 
% figure (5)
% surf(W)

start2 = 2000;
start2 = 0;
figure(1);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), mean_Vx(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x)');xlabel('t, ps');ylabel('m/s');
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), int_Vx_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x) dt');xlabel('t, ps');ylabel('nm');

% sumFxup = cumsum(M(:,cFx_wk_up));
% sumFxdw = cumsum(M(:,cFx_wk_dw));
% 
% len = size(sumFxup,1);
% N = 10
% wind = fix(len / N);
% wind_05 = fix(wind / 2);
% 
% ma_sumFxup = zeros(size(sumFxup));
% ma_sumFxdw = zeros(size(sumFxdw));
% 
% for i = 1 : len - wind
%     ma_sumFxup(i + wind_05) = mean(sumFxup(i:i+wind));
%     ma_sumFxdw(i + wind_05) = mean(sumFxdw(i:i+wind));
% end
% 
% 
% figure(2);
% subplot(2,1,1);
% plot (0.001*M(:,1), ma_sumFxup, [0.001*M(1,1) 0.001*M(end,1)],[0 0],'black');title('sum F_x _w_k up'); xlabel('t, ps');
% subplot(2,1,2);
% plot (0.001*M(:,1), ma_sumFxdw, [0.001*M(1,1) 0.001*M(end,1)],[0 0],'black');title('sum F_x _w_k down'); xlabel('t, ps');
%  
start1 = 40000;
start1 = 1;
figure(3);
subplot(2,1,1);
plot (0.001*M(start1:end,1), cumsum(M(start1:end,cFx_wk_up)), [0.001*M(1,1) 0.001*M(end,1)],[0 0],'black');title('cumsum F_x _w_k up'); xlabel('t, ps');
subplot(2,1,2);
plot (0.001*M(start1:end,1), cumsum(M(start1:end,cFx_wk_dw)), [0.001*M(1,1) 0.001*M(end,1)],[0 0],'black');title('cumsum F_x _w_k down'); xlabel('t, ps');
 

mean_f_wk_up = cumsum(M(:,cFx_wk_up))./((1:size(M,1))');
mean_f_wk_down = cumsum(M(:,cFx_wk_dw))./((1:size(M,1))');
start2 = 5000;

% start2 = 0;

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
 












start1 = 31000;
start1 = 1;
mean_dF_wk_x = cumsum(M(start1:end,cFx_wk_up)-M(start1:end,cFx_wk_dw))./((1:size(M,1)-start1+1)');

figure(6);
subplot(3,1,1);
start2 = 1000;
% start2 = 0;
plot (0.001*M(start1+start2:end,1), mean_dF_wk_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0],'black');title('mean (F_x _w_k up - F_x _w_k down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,2);
start2 = 10000;
plot (0.001*M(start1+start2:end,1), mean_dF_wk_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0],'black');title('mean (F_x _w_k up - F_x _w_k down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,3);
start2 = 20000;
plot (0.001*M(start1+start2:end,1), mean_dF_wk_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0],'black');title('mean (F_x _w_k up - F_x _w_k down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
% 


start1 = 1;
start2 = 2000;
start2 = 0;

figure(7);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), 1000*M(start1+start2:end,cVx_wk_up), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('V_x _w_k up');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), 1000*M(start1+start2:end,cVx_wk_dw), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('V_x _w_k down');xlabel('t, ps');ylabel('m/s');
 

start1 = 1;
% start1 = 100000;
start2 = 50000;
start2 = 5000;

mean_V_wk_x_up = 1000*cumsum(2*(M(start1:end,cVx_wk_up))/natoms_worked)./((1:size(M,1)-start1+1)');
int_V_wk_x_up_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_wk_up))/natoms_worked);



figure(8);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), mean_V_wk_x_up(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x _w_k up)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), mean_V_wk_x_up(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x _w_k up) dt');xlabel('t, ps');ylabel('nm');
 

mean_V_wk_dw_x = 1000*cumsum(2*(M(start1:end,cVx_wk_dw))/natoms_worked)./((1:size(M,1)-start1+1)');
int_V_wk_x_dw_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_wk_dw))/natoms_worked);



figure(9);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), mean_V_wk_dw_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x _w_k down)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), int_V_wk_x_dw_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x _w_k down) dt');xlabel('t, ps');ylabel('nm');
 


mean_dV_wk_x = 1000*cumsum(2*(M(start1:end,cVx_wk_up)-M(start1:end,cVx_wk_dw))/natoms_worked)./((1:size(M,1)-start1+1)');
int_dV_wk_x_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_wk_up)-M(start1:end,cVx_wk_dw))/natoms_worked);



figure(10);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), mean_dV_wk_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x _w_k up - V_x _w_k down)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), int_dV_wk_x_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x _w_k up - V_x _w_k down) dt');xlabel('t, ps');ylabel('nm');
 

% start1 = 1;
% start2 = 50000;
% start2 = 0;

mean_V_x_up = 1000*cumsum(2*(M(start1:end,cVx_up))/natoms)./((1:size(M,1)-start1+1)');
int_V_x_up_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_up))/natoms);



figure(12);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), mean_V_x_up(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x up)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), mean_V_x_up(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x up) dt');xlabel('t, ps');ylabel('nm');
 

mean_V_dw_x = 1000*cumsum(2*(M(start1:end,cVx_dw))/natoms)./((1:size(M,1)-start1+1)');
int_V_x_dw_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_dw))/natoms);



figure(13);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), mean_V_dw_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x down)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), int_V_x_dw_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x down) dt');xlabel('t, ps');ylabel('nm');
 


mean_dV_x = 1000*cumsum(2*(M(start1:end,cVx_up)-M(start1:end,cVx_dw))/natoms)./((1:size(M,1)-start1+1)');
int_dV_x_dt = 1000*tstep*1e9*cumsum(2*(M(start1:end,cVx_up)-M(start1:end,cVx_dw))/natoms);



figure(14);
subplot(2,1,1);
plot (0.001*M(start1+start2:end,1), mean_dV_x(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (V_x up - V_x down)');xlabel('t, ps');ylabel('m/s');
 
subplot(2,1,2);
plot (0.001*M(start1+start2:end,1), int_dV_x_dt(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('int (V_x up - V_x down) dt');xlabel('t, ps');ylabel('nm');
 


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

% figure(11);
% plot (1 : N, means);



start1 = 1047;
dVx = 1000*2*(M(start1:end,cVx_wk_up)-M(start1:end,cVx_wk_dw))/natoms_worked;

len = size(dVx,1);
means = [];
N = 20000
wind = fix(len / N);
for i = 1 : N
    m = mean(dVx(1 + (i-1)*wind:i*wind));
    means = [means m];
end
mat_ozh = mean(means)
D = sum((means - mat_ozh).^2)/(size(means, 2)-1);
sigma = sqrt(D / size(means, 2))

