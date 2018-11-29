addpath('C:\Users\lamar\Documents\ECA UFSC\Projeto Daniel\FalcOpt');

A = [0.7887 0       0.2820  0
     0      0.8486  0       0.2155
     0      0       0.6808  0 
     0      0       0       0.7654];
 
 B = [  0.6439  0.1823
        0.0973  0.4836
        0       0.9657
        0.7362  0       ];
  
Nx = 4;
Nu = 2;
N  = 15;
Ts = 15;

variable_stepSize.active = false;

J.Q = [ eye(Nx,Nx)  ,zeros(Nx,Nu); 
        zeros(Nu,Nx),zeros(Nu,Nu)];
J.R = eye(Nu);
J.P = J.Q;
J.trackReference = true;

dumin = [-0.5;-0.5];
dumax = [ 0.5; 0.5];

eps             = 1e-3;      % tolerance
merit_function  = 0;         % merit function
debug           = 3;        % level of debug

gradients   = 'matlab';
precision   = 'double';
contractive = false;        % no constractive constraints
terminal    = false;        % no terminal constraints

A_ = [A,B;zeros(Nu,Nx),eye(Nu)];
B_ = [B;eye(Nu)];

dynamics = @(x,du) A_*x+B_*du;

constraint = {};

nn = 0; % Número de restrições não-lineares

info = falcopt.generateCode(dynamics, N, Nx+Nu, Nu, J              ,...
                'variable_stepSize' , variable_stepSize         ,...
                'constraints_handle', constraint                ,...
                'nn'                , nn                        ,...
                'gradients'         , gradients                 ,...
                'lb'    , dumin                      ,...
                'ub'    , dumax                      ,...
                'contractive'       , contractive               ,...
                'terminal'          , terminal                  ,...
                'debug'             , debug                     ,...
                'merit_function'    , merit_function            ,...
                'eps'               , eps                       ,...
                'precision'         , precision                 ,...
                'name'              , 'controladorFalcOpt'   ,...
                'gendir'            , 'generatedCode'           ,...
                'buildTypes'        , {'mex', 'production'},...[]
                'maxIt'             , 8000                 ,...
                'maxItLs'           , 50);
