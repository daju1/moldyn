clear;
close all;

[FILENAME, PATHNAME] = uigetfile(['.txt']);
cd(PATHNAME);

M = dlmread([PATHNAME, FILENAME], ',');

mean(M,1)

figure(1);
subplot(2,1,1);
plot (0.001*M(:,1), M(:,4));title('sum V_x'); xlabel('t, ps');
subplot(2,1,2);
plot (0.001*M(:,1), cumsum(M(:,4)));title('cumsum V_x'); xlabel('t, ps');
 
figure(2);
subplot(2,1,1);
plot (0.001*M(:,1), M(:,7));title('sum F_x up'); xlabel('t, ps');
subplot(2,1,2);
plot (0.001*M(:,1), M(:,10));title('sum F_x down'); xlabel('t, ps');
 
figure(3);
subplot(2,1,1);
plot (0.001*M(:,1), cumsum(M(:,7)));title('cumsum F_x up'); xlabel('t, ps');
subplot(2,1,2);
plot (0.001*M(:,1), cumsum(M(:,10)));title('cumsum F_x down'); xlabel('t, ps');
 

mean_fup = cumsum(M(:,7))./((1:size(M,1))');
mean_fdown = cumsum(M(:,10))./((1:size(M,1))');
start2 = 50000;

figure(4);
subplot(2,1,1);
plot (0.001*M(1+start2:end,1), mean_fup(1+start2:end), [0.001*M(1+start2,1) 0.001*M(end,1)],[0 0]);title('mean F_x up'); xlabel('t, ps');ylabel('10^1^2 N/mol');
subplot(2,1,2);
plot (0.001*M(1+start2:end,1), mean_fdown(1+start2:end), [0.001*M(1+start2,1) 0.001*M(end,1)],[0 0]);title('mean F_x down'); xlabel('t, ps');ylabel('10^1^2 N/mol');
 










start1 = 2500;
mean_t = cumsum(M(start1:end,7)-M(start1:end,10))./((1:size(M,1)-start1+1)');

figure(5);
subplot(3,1,1);
start2 = 1000;
plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,2);
start2 = 20000;
plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,3);
start2 = 50000;
plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 












start1 = 1;
mean_t = cumsum(M(start1:end,7)-M(start1:end,10))./((1:size(M,1)-start1+1)');

figure(6);
subplot(3,1,1);
start2 = 1000;
plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,2);
start2 = 20000;
plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 
subplot(3,1,3);
start2 = 70000;
plot (0.001*M(start1+start2:end,1), mean_t(1+start2:end), [0.001*M(start1+start2,1) 0.001*M(end,1)],[0 0]);title('mean (F_x up - F_x down)');xlabel('t, ps');ylabel('10^1^2 N/mol');
 

start1 = 2500;
dF = M(start1:end,7)-M(start1:end,10);
dF = M(start1:end-10000,7)-M(start1:end-10000,10);
len = size(dF,1);
means = [];
N = 20
wind = fix(len / N);
for i = 1 : N
    m = mean(dF(1 + (i-1)*wind:i*wind));
    means = [means m];
end
mat_ozh = mean(means)
D = sum((means - mat_ozh).^2)/(size(means, 2)-1);
sigma = sqrt(D / size(means, 2))

