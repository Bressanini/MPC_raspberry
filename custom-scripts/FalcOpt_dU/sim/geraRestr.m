
tipo = 'lin';

NX = 4;
NU = 2;

NIT = 100;

Q = [1.0,1.0,1.0,1.0];
P = Q;
R = [1.0,1.0];

u0 = [  3.15;
        3.15];
    
x0 = [  12.6;
        13.0;
        4.80;
        4.90];
x0sim = [   x0(1) + 5;
            x0(2) - 5;
            x0(3);
            x0(4)];

uMIN = u0/10;
uMAX = u0*2;
    
xMIN = x0/10;
    
xMAX = x0*2;

uREF = [u0*ones(1,40), [4.0;2.5]*ones(1,60)];

xREF = [x0*ones(1,40), [12.0053;15.2714;2.8335;7.5676]*ones(1,60)];
    
%% Fim da entrada de dados

if tipo == 'lin'
        
    uMIN = uMIN - u0;
    uMAX = uMAX - u0;
    xMIN = xMIN - x0;
    xMAX = xMAX - x0;
    uREF = uREF - u0*ones(1,NIT);
    xREF = xREF - x0*ones(1,NIT);
    x0sim = x0sim - x0;

end

csvwrite('Q.csv'    ,Q);
csvwrite('P.csv'    ,P);
csvwrite('R.csv'    ,R);
csvwrite('uMIN.csv' ,uMIN);
csvwrite('uMAX.csv' ,uMAX);
csvwrite('xMIN.csv' ,xMIN);
csvwrite('xMAX.csv' ,xMAX);
csvwrite('uREF.csv' ,uREF);
csvwrite('xREF.csv' ,xREF);
csvwrite('u0.csv'   ,u0);
csvwrite('x0.csv'   ,x0);
csvwrite('x0sim.csv',x0sim);

