% Programa para plotar resultados

u = csvread('U.csv'); % Lendo as manipuladas
u = u(:,1:end-1); % Por colocar vírgulas no final do csv, ele entende que existe uma última coluna com zeros

x = csvread('X.csv');
x = x(:,1:end-1);

uREF = csvread('uREF.csv');

xREF = csvread('xREF.csv');


figure;

ax(1) = subplot(2,1,1);
stairs(1:(length(x)-1),xREF','--');
hold on;
stairs(1:length(x),x');
xlabel('Iteração'); ylabel('Nível / cm')
legend('h_{1,ref}','h_{2,ref}','h_{3,ref}','h_{4,ref}',...
    'h_1','h_2','h_3','h_4',...
    'Location','bestoutside')

ax(2) = subplot(2,1,2);
stairs(1:(length(x)-1),uREF','--');
hold on;
stairs(2:length(x),u');
xlabel('Iteração'); ylabel('Tensão / V')
legend('v_{1,ref}','v_{2,ref}',...
    'v_1','v_2',...
    'Location','bestoutside')

linkaxes(ax,'x')
axis([1 length(x) -inf inf])