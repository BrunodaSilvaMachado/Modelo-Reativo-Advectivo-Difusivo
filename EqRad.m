% Copyright (C) 2023 B da Silva Machado, GB Alvarez, DC Lobão
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%==========================================================================
%
%       Reactive-Advective-Diffusive Models for the Growth of Gliomas 
%                        Treated with Radiotherapy
%                B da Silva Machado, GB Alvarez, DC Lobão
%        Semina: Ciências Exatas e Tecnológicas, 2023, v44, e47321-e47321
%           DOI: https://doi.org/10.5433/1679-0375.2023.v44.47321
%
% Version 000:
% Implemented by Bruno Machado da Silva and Gustavo Benitez Alvarez
% UFF - Volta Redonda, RJ, Brazil
% Oct 19th, 2023
%
%==========================================================================
function [] = EqRad(frac_i, varargin)
%EQRAD Compara as equações de crescimento tumoral diante a radioterapia de rockne e steins.
% Parâmetro obrigatório: 
%   frac_i = <valor>; Esquema de fracionamento com níveis de 1 até 8.
% Parâmetros opcionais:
%   input = {D, <valor>, rho, <valor>, V, <valor>, Lx,<valor>, alpha, <valor>, beta, <valor>};
%   Parâmetros do modelo de crescimento tumoral.
%   days = 80; Número de days na simulação.
%   position = 0; posição inicial do tumor cerebral.
%   nos = 1400; Número de nós na malha
%   cinitial = 1; Seleciona a condição inicial. Atualmente disponível as
%   condições inicias 1 ou 2.
%   directory = <CAMINHO>; Diretório para save os arquivos de saída
%   models = {'rd-linear','rd','rad'}; Modelos para crescimento tumoral com
%   radioterapia. Exemplos e descrições abaixo.
%   graphs = {'ci','ci-comp','crm','crm-comp','none'}; Opções de gráficos
%   de saída. Exemplos e descrições abaixo.
%   view = <'on', 'off'>; Se 'on' os gráficos gerados são exibidos na tela.
%   save = <true, false>; Se true salva a imagem no diretório especificado.
%   wt = <true, false>; Se verdadeiro grava uma tabela com os resultados em um arquivo txt. 
%   ext = '-dpng'; Extensão do ficheiro de saída. Usado quando save ou wt são verdadeiros.
%

%% Exemplos 
% Entrada:
% ab = 10; % sensibilidade do  tecido ao fracionamento da dose
% alpha = 0.0305;
% input = { 'D', 3.9e-5/24,...     % mobilidade celular cm^2/dia convertido em horas 
%           'rho', 0.0453/24,...   % taxa de proliferação 1/dia convertido em horas
%           'V', 0.01/24,...       %cm/dia convertidos em horas
%           'Lx', 20,...           % comprimento do domínio (cérebro) cm
%           'alpha', alpha,...     % parâmetro de sensibilidade à radiação Gy^-1
%           'beta', alpha/ab
% }
% 
%
% Modelos:
% models = {'rd-linear',... % Modelo de crescimento tumoral de Rockne e Swason
%                           % Tipo Reativo-difusivo
%           'rd',...        % Modelo de Rockne, mas crescimento logístico. Tipo
%                           % Reativo-difusivo.
%           'rad'           % Nosso modelo. Baseado no modelo de Steins. Tipo
%                           % Reativo-Adivectivo-difusivo.
%           }; 
%
%
% Gráficos:
% Os gráficos impressos por este programa são: 
% graphs = {
%   'ci',...       % condição inicial
%   'ci-comp',...  % comparação entre as condições iniciais dos modelos 'rd',
%                  % 'rd-linar', 'rad'. Deve ter pelo menos 2 modelos selecionados para usar esta opção.
%   'cl',...       % concentração final de células tumorais
%   'cl-comp',...  % comparação entre as concentrações finais de células tumorais dos modelos 'rd',
%                  % 'rd-linar', 'rad'. Deve ter pelo menos 2 modelos selecionados para usar esta opção.
%   'crm',...      % concentração final, ratio e concentração máxima das células tumorais
%   'crm-comp',... % comparação entre as concentrações finais, ratios e concentrações máximas das céluas
%                  % tumorais dos modelos 'rd', 'rd-linar', 'rad'. Deve ter pelo menos 2 modelos selecionados para usar esta opção.
%   'none'         % Não imprimir gráficos
%

%EQRAD Compares tumor growth equations to rockne and steins radiotherapy.
% Required parameter: 
%   frac_i = <value>; Fractionation scheme with levels from 1 to 8.
% Optional parameters:
%   input = {D, <value>, rho, <value>, V, <value>, Lx,<value>, alpha, <value>, beta, <value>};
%   Parameters of the tumor growth model.
%   days = 80; Number of days in the simulation.
%   position = 0; Initial position of the brain tumor.
%   nos = 1400; Number of nodes in the mesh
%   cinitial = 1; Selects the initial condition. Currently available
%   initial conditions 1 or 2.
%   directory = <PATH>; Directory for saving output files
%   models = {'rd-linear','rd','rad'}; Models for tumor growth with
%   radiotherapy. Examples and descriptions below.
%   graphs = {'ci','ci-comp','crm','crm-comp','none'}; Graph options
%   output. Examples and descriptions below.
%   view = <'on', 'off'>; If 'on' the generated graphics are displayed on the screen.
%   save = <true, false>; If true saves the image in the specified directory.
%   wt = <true, false>; If true saves a table with the results in a txt file. 
%   ext = '-dpng'; Extension of the output file. Used when save or wt are true.
%

%% Examples 
% Input:
% ab = 10; % tissue sensitivity to dose fractionation
% alpha = 0.0305;
% input = { 'D', 3.9e-5/24,...     % cell mobility cm^2/day converted to hours 
%           'rho', 0.0453/24,...   % proliferation rate 1/day converted to hours
%           'V', 0.01/24,...       %cm/day converted to hours
%           'Lx', 20,...           % length of domain (brain) cm
%           'alpha', alpha,...     % radiation sensitivity parameter Gy^-1
%           'beta', alpha/ab
% }
% 
%
% Models:
% models = {'rd-linear',... % Rockne and Swason tumor growth model
%                           % Reactive-diffusive type
%           'rd',...        % Rockne model, but logistic growth. Type
%                           % Reactive-diffusive.
%           'rad'           % Our model. Based on the Steins model. Type
%                           % Reactive-Advective-diffusive.
% }; 
%
%
%
% Graphics:
% The graphs printed by this program are: 
% graphs = {
%               'ci',...       % initial condition
%               'ci-comp',...  % comparison between the initial conditions of the 'rd' models,
%                              % 'rd-linar', 'rad'. You must have at least 2 models selected to use this option.
%               'cl',...       % final concentration of tumor cells
%               'cl-comp',...  % comparison between the final tumor cell concentrations of the 'rd' models,
%                              % 'rd-linar', 'rad'. You must have at least 2 models selected to use this option.
%               'crm',...      % final concentration, radius and maximum concentration of tumor cells
%               'crm-comp',... % comparison of final concentrations, radii and maximum concentrations of tumor cells
%                              % tumor cells of the 'rd', 'rd-linar', 'rad' models. 
%                              % You must have at least 2 models selected to use this option.
%               'none'         % Do not print graphs
%           }
%

ent = opt(frac_i, varargin{:});
D = ent.inputStruct.D/24;
rho = ent.inputStruct.rho/24;
V = ent.inputStruct.V/24;
Lx = ent.inputStruct.Lx;
alpha = ent.inputStruct.alpha;
beta = ent.inputStruct.beta;
%==========================================================
%% Discretize the spatial and temporal domain

N = ent.nos; % Number of nodes in the mesh
ci0 = ent.cinitial;
N_1 = N + 1;
delx = Lx/N; % Spatial variation (discretization)
delx_a = delx/Lx;
x = (0:delx:Lx)';
Xa = x/Lx;
Xc = ent.position/Lx;
if (Xc > 0.5) 
    V = -V; 
end
Nt = ent.days*24; %Number of time steps. 1920 hours = 80 days
delt = 1; % Temporal variation
delt_a = delt * 0.0011; %delt_a = delt * rho; % Delta T adimensionalizado
Limiar = 0.6126*Lx^(3); %Limiar is used to calculate the radius
U_0 = Limiar*exp(100*(1.41/Lx)^2); %U_0 which generates a radius of 1.41 cm 
U_max = 4.2e8; %maximum cell concentration/cm^3 (Steins).
%% Normalização
lambda = 0.5*(D/(rho*Lx^2))*(delt_a/delx_a^2); %Lambda
vi = V/(rho *Lx) *(delt_a/delx_a); % advective term
l = lambda * ones(N_1,1);
vi_p = max(vi,0);
vi_m = min(vi,0);
s_p = l + vi_p;
s_m = l - vi_m;
tau = delt_a/2;% Constant of the non-linear term
%% ================================================
%% Dose fractionation
tratament = [
    struct('dose',(0) , 'fraction',0 , 'reforce',false ); %without treatment
    struct('dose',[60,4.2] , 'fraction',1 , 'reforce',true );
    struct('dose',(20), 'fraction',3 , 'reforce',false );
    struct('dose',[12.2,3.2] , 'fraction',5 , 'reforce',true );
    struct('dose',[6.0,6.0], 'fraction',10 , 'reforce',false );
    struct('dose',[2.8,3.5,6.5] , 'fraction',15 , 'reforce',false );
    struct('dose',[2.0,2.9,2.9,3.0,2.0] , 'fraction',25 , 'reforce',false );
    struct('dose',[1.8,1.8,1.8,1.8,1.8,1.8,1.8] , 'fraction',35 , 'reforce',false )
    ];
trat = tratament(frac_i);
tap = fractionation(trat.fraction,trat.reforce);
dailydose = dosage(trat.dose,trat.fraction,trat.reforce);
%% ================================================
%% Initial condition
initialCondition = [
    initialConditionExp(Lx, U_max, Xa, Xc),...
    initialConditionC1(U_0, Limiar, 0,0.6/Lx,1.41/Lx,2.0/Lx, N_1, U_max, Xa, Xc)
    ];
R_0 = sum(U_max * initialCondition(:,ci0) > Limiar) * delx;

% Models
rd = @(d,t,H_n) rockne(alpha, beta, rho, d, t,...
        U_max, delx, Limiar, tau, l, H_n, N_1);
rdn = @(d,t,U_n) rockne_n_linear(alpha, beta, rho, d, t,...
        U_max, delx, Limiar, tau, l, U_n, N_1);
rad = @(d,t,W_n) stein_modificado(alpha, beta, rho, d, t,...
        U_max, delx, Limiar, tau,l,s_p,s_m,W_n,N_1);
    
models = {};
% Reative-Difusive (linear)
if any(strcmp('rd-linear', ent.models))
    H_n = initialCondition(:,ci0);
    models{end+1} = struct('edf',rd,'A_n',H_n,'A_i',U_max*H_n,'concentration',...
        zeros(1,Nt),'ratio',zeros(1,Nt),'tagTitle','RD linear','tagFile','rdLinear');
end
%Reative-Difusive 
if any(strcmp('rd', ent.models))
    U_n = initialCondition(:,ci0);
    models{end+1} = struct('edf',rdn,'A_n',U_n,'A_i',U_max*U_n,'concentration',...
        zeros(1,Nt),'ratio',zeros(1,Nt),'tagTitle','RD','tagFile','rd');
end
% Reative-Advective-Difusive
if any(strcmp('rad', ent.models))
    W_n = initialCondition(:,ci0);
    models{end+1} = struct('edf',rad,'A_n',W_n,'A_i',U_max*W_n,'concentration',...
        zeros(1,Nt),'ratio',zeros(1,Nt),'tagTitle','RAD','tagFile','rad');
end
tic;
%figure(1);
for n = 1:Nt %n =1 is the Initial Condition
    %Therapy Application
    isTap = (tap == n);
    d = dailydose(find(isTap,1));
    therapy = sum(isTap);
    %Solving the dimensionless system
    for i = 1:length(models)
        [A_n_1,cn, rn] = models{i}.edf(d, therapy, models{i}.A_n);
        models{i}.A_n = A_n_1;
        models{i}.concentration(n) = cn;
        models{i}.ratio(n) = rn;
    end
end
%End of the time loop
%% ================================================
toc
% Graphics
% initial concentration
if any(strcmp('ci', ent.graphics))
    for i = 1:length(models)
       graphicsCI(ent,trat.fraction,x, models{i}.A_i,models{i}.tagTitle,...
           models{i}.tagFile);
    end
end
% initial concentration comparison
if any(strcmp('ci-comp', ent.graphics))
    n = length(models);
    m = length(models{1}.A_i);
    A_i = zeros(m,n);
    ll = {};
    for i = 1:n
        A_i(:,i) = models{i}.A_i;
        ll{i} = models{i}.tagTitle;
    end
    graphicsCIComp(ent,trat.fraction,x, A_i,ll);
end
% final concentration
if any(strcmp('cl', ent.graphics))
    for i = 1:length(models)
    graphicsCL(ent, trat.fraction, models{i}.tagTitle, models{i}.tagFile,...
        U_max,x, models{i}.A_n, Limiar);
    end
end
% final concentration comparison
if any(strcmp('cl-comp', ent.graphics))
    n = length(models);
    m = length(models{1}.A_n);
    A_n = zeros(m,n);
    ll = {};
    for i = 1:n
        A_n(:,i) = U_max*models{i}.A_n;
        ll{i} = models{i}.tagTitle;
    end
    graphicsCLComp(ent,trat.fraction,x, A_n,ll);
end
% concentration, radius and maximum concentration
if any(strcmp('crm', ent.graphics))
    for i = 1:length(models)
    graphicsCRM(ent, trat.fraction, models{i}.tagTitle, models{i}.tagFile,...
        U_max,x, models{i}.A_i, models{i}.A_n, R_0, models{i}.ratio,...
        models{i}.concentration, Limiar);
    end
end
% concentration, radius and maximum concentration comparison
if any(strcmp('crm-comp', ent.graphics))
    n = length(models);
    m_c = length(models{1}.concentration);
    m_r = length(models{1}.ratio);
    c_i = zeros(m_c,n);
    r_i = zeros(m_r,n);
    ll = {};
    for i = 1:n
        c_i(:,i) = models{i}.concentration;
        r_i(:,i) = models{i}.ratio;
        ll{i} = models{i}.tagTitle;
    end
    graphicsRMComp(ent,trat.fraction, ll, U_max, R_0, r_i, c_i);
end

if ent.wt
    n = length(models);
    m_c = length(models{1}.concentration);
    m_r = length(models{1}.ratio);
    m_a = length(models{1}.A_n);
    c_i = zeros(m_c,n);
    r_i = zeros(m_r,n);
    A_n = zeros(m_a,n);
    for i = 1:n
        c_i(:,i) = U_max*models{i}.concentration;
        r_i(:,i) = models{i}.ratio;
        A_n(:,i) = U_max*models{i}.A_n;
    end
    
    writeStructModel(ent,trat.fraction, x, R_0, r_i, c_i, A_n)
end
end

function [results] = opt(fraction, varargin)
    cellstr = @(x) any(cellfun(@ischar, x));
    defaultDays = 80;
    defaultPosition = 0;
    defaultExt = '-dpng';
    defaultNos = 1400;
    defaultCInitial = 1;
    defaultDirectory = './';
    defaultSave = true;
    defaultWriteTable = false;
    defaultModels = {'rd','rd-linear','rad'};
    defaultGraphics = {'ci','ci-comp','cl','crm','crm-comp'};
    defaultView = 'on';
    defaultInput = {'D',3.9e-5,'rho',0.0453,'V',0.01,'Lx',20,...
        'alpha',0.0305,'beta',0.00305};
    
    p = inputParser;
    p.KeepUnmatched = true;
    addRequired(p, 'fraction',@isnumeric);
    addParameter(p, 'days',defaultDays, @isnumeric);
    addParameter(p, 'position', defaultPosition, @isnumeric);
    addParameter(p, 'ext', defaultExt,@ischar);
    addParameter(p, 'save', defaultSave,@islogical);
    addParameter(p, 'nos', defaultNos, @isnumeric);
    addParameter(p, 'cinitial',defaultCInitial , @isnumeric);
    addParameter(p, 'directory',defaultDirectory,@ischar);
    addParameter(p, 'wt', defaultWriteTable,@islogical);
    addParameter(p, 'models', defaultModels,cellstr);
    addParameter(p, 'graphics', defaultGraphics,cellstr);
    addParameter(p, 'view', defaultView,@ischar);
    addParameter(p, 'input', defaultInput,@iscell);
    parse(p, fraction, varargin{:});
    
    e = inputParser;
    e.KeepUnmatched = true;
    addParameter(e, 'D',defaultInput{2}, @isnumeric);
    addParameter(e, 'rho',defaultInput{4}, @isnumeric);
    addParameter(e, 'V',defaultInput{6}, @isnumeric);
    addParameter(e, 'Lx',defaultInput{8}, @isnumeric);
    addParameter(e, 'alpha',defaultInput{10}, @isnumeric);
    addParameter(e, 'beta',defaultInput{12}, @isnumeric);
    parse(e,p.Results.input{:});
    
    results = p.Results;
    results.inputStruct = e.Results;
end

function [ci] = initialConditionExp(Lx, U_max, Xa, Xc)
    ci = (Lx^3/U_max)*exp(-100*(Xa - Xc).^2); %normalized initial condition;
end
function [ci] = initialConditionC1(C_m ,C_l, C_z, Xm, Xl, Xz, N_1, M, Xa, Xc)
    if(Xc ~= 0)
        disp('Xc diferente de zero. Sem suporte a deslocamento');
    end
    Xac = Xa;
    ci = zeros(N_1,1);
    id = (Xac <= Xm);
    ci(id) = C_m;
    id = (Xac > Xm & Xac <= Xl);
    ci(id) = ((C_m - C_l)/(Xm - Xl))*(Xac(id) - Xm) + C_m;
    id = (Xac > Xl & Xac <= Xz);
    ci(id) = ((C_l - C_z)/(Xl - Xz))*(Xac(id) - Xl) + C_l;
    id = Xac > Xz;
    ci(id) = 0;
    ci = ci/M;
end

function [H_n_1, h, r] = rockne(alpha, beta, rho, dose, therapy,...
    MU,delx, Limiar, tau, l, H_n, N_1)
% Reative difusive
g = tau * XRT(alpha, beta, rho, dose, therapy);
%Calculation of p and q Dimensionaless
p = 1 + 2*l - g;
q = 1 - 2*l + g;
%End of phi and psi calculation
%Assembling the MA and ME matrices
MA = spdiags([-l p -l], [-1 0 1], N_1, N_1);
MA(1,1:3) =  [-3,4,-1];
MA(N_1,N_1-2:N_1) = [1,-4,3];
ME = spdiags([l q l], [-1 0 1], N_1, N_1);
ME(1,1:3) = [3,-4,1];
ME(N_1,N_1-2:N_1) = [-1,4,-3];
%Solving system A dimensionless
b = ME*H_n;
H_n_1 = MA\b;
%maximum concentration
h = max(H_n_1);
%Calculating the tumor radius without dimensionless
%X_c = 0 is the radius, not the diameter.
tam = sum(MU*H_n_1 > Limiar);
r = 0;
if tam > 0
    r = tam * delx;
end
end

function [U_n_1, u, r] = rockne_n_linear(alpha, beta, rho, dose, therapy,...
    MU,delx, Limiar,tau,l, U_n, N_1)
% Non-linear Reative difusive
g = tau * XRT(alpha, beta, rho, dose, therapy);
% Calculation of p and q Dimensionaless
p = 1 + 2*l - g + 2 * tau * U_n;
q = 1 - 2*l + g;
%MA
MA = spdiags([-l p -l], [-1 0 1], N_1, N_1);
MA(1,1:3) =  [-3,4,-1];
MA(N_1,N_1-2:N_1) = [1,-4,3];
%MD
MD = spdiags([l q l], [-1 0 1], N_1, N_1);
MD(1,1:3) =  [3,-4,1];
MD(N_1,N_1-2:N_1) = [-1,4,-3];
%RD
b_u = MD*U_n;
U_n_1 = MA\b_u;
u = max(U_n_1);
%Calculating the tumor radius without dimensionless
%X_c = 0 is the radius, not the diameter.
tam = sum(MU*U_n_1 > Limiar);
r = 0;
if tam > 0
    r = tam * delx;
end
end

function [W_n_1,wm,r] = stein_modificado(alpha, beta, rho, dose, therapy,...
    MU, delx,Limiar,tau,l, s_p, s_m, W_n, N_1)
% Reative Advctive Difusive
g = tau * XRT(alpha, beta, rho, dose, therapy);
p_w = 1 + 2*l - g + 2 * tau * W_n;
r = 1 + g - s_p - s_m;
%MA
MA = spdiags([-l p_w -l], [-1 0 1], N_1, N_1);
MA(1,1:3) =  [-3,4,-1];
MA(N_1,N_1-2:N_1) = [1,-4,3];
%ME
ME = spdiags([s_p r s_m], [-1 0 1], N_1, N_1);
ME(1,1:3) =  [3,-4,1];
ME(N_1,N_1-2:N_1) = [-1,4,-3];
%RAD
b_w = ME*W_n;
W_n_1 = MA\b_w;
wm = max(W_n_1);
%Calculating the tumor radius without dimensionless
tam = sum(MU*W_n_1 > Limiar);
r = 0;
if tam > 0
    r = tam * delx;
end
end

function [ret] = XRT(alpha, beta, rho, dose, therapy)
% Radiotherapy equal doses without boost, delta t constant
ret = 1;

if (therapy)
    % Probability S of cell survival with the effective biological dose
    BED = alpha * dose + beta * dose^2;
    S = exp(-BED);
    ret = (1 - (1 - S)/rho);
end
end


function [tap] = fractionation(fraction, reforce)
% Treatment equal dose in Gy with booster
N = fraction + 1*(reforce == true);
tap = zeros(1, N); % application time in hours
tap(1) = 8;% hours
DAYS = 5; %treatment days
if (N == 0) % No Treatment
    tap(1)  = 0;
else
    for i = 1: N - 1
        % for tap longer than a week, add a 72-hour break
        % between applications.
        tap(i + 1) = tap(i) + (mod(i,DAYS) ~= 0) * 24 + (mod(i,DAYS) == 0) * 72;
    end
end
end

function [ndoses] = dosage(doseweek, fid, reforce)
    % Transforms weekly dose into daily dose or with or without booster
    N = fid + 1*(reforce == true);
    DAYS = 5; %treatment days
    index = 1;
    ndoses = zeros(1, N);
    if (fid == 0) % No dosage
        ndoses(1) = 0;
    else
        for i = 1:fid
            ndoses(i) = doseweek(index);
            if(mod(i,DAYS) == 0)
                index = index + 1;
            end
        end
    end
    
    if(reforce)
        ndoses(N) = doseweek(end);
    end
end

function writeStructModel(ent, fraction ,x, R_0, r_i, c_i, A_n)
% Prints graphs of tumor concentration, radius and maximum concentration.
%
DIRECTORY = ent.directory; % Output directory
TYPE_TRATAMENT = sprintf('Dot%02ddays%dX%dCi%d', ...
        fraction, ent.days,ent.position, ent.cinitial);
ext = ent.ext;
Nt = ent.days * 24;
eixoT = linspace(0, ent.days, Nt);
%final concentration
writetable(array2table([x,A_n]),strcat(DIRECTORY,'mConcentracao',TYPE_TRATAMENT,'.',ext));

% evolution of the glioma radius
writetable(array2table([eixoT',r_i,R_0*ones(Nt,1)]),strcat(DIRECTORY,'mRaioGlioma',TYPE_TRATAMENT,'.',ext));

%maximum concentration
%lc_i =  log(c_i);
writetable(array2table([eixoT',c_i]),strcat(DIRECTORY,'mConcentracaoMax',TYPE_TRATAMENT,'.',ext));

end

function graphicsCIComp(ent,fraction, x, AB, clegend)
%Initial Condition Figure
%
fig1 = figure(1); set(fig1,'visible',ent.view);
clf
hold on
for A = AB
    plot (x, A,'LineWidth',1.5);
end
hold off
legend(clegend);
set(gca,'FontSize',12);
ylabel('Concentração (células/cm)');
xlabel('X (cm)');
title('Condição inicial');
if ent.save
    DIRECTORY = ent.directory;% Output directory
    TYPE_TRATAMENT = sprintf('Dot%02ddays%dX%dCi%d', ...
        fraction, ent.days,ent.position, ent.cinitial);
    ext = ent.ext;
    saveas (gcf, strcat(DIRECTORY,'rdcConcentracaoInicial',TYPE_TRATAMENT), ext);
end
end

function graphicsCI(ent,fraction, x, AB, tagTitle, tagFile)
%Initial Condition Figure
%
fig1 = figure(); set(fig1,'visible',ent.view);
clf
plot (x, AB,'LineWidth',1.5);
set(gca,'FontSize',12);
ylabel('Concentração (células/cm)');
xlabel('X (cm)');
title(sprintf('Condição inicial (%s)', tagTitle));
if ent.save
    DIRECTORY = ent.directory;% Output directory
    TYPE_TRATAMENT = sprintf('Dot%02ddays%dX%dCi%d', ...
        fraction, ent.days,ent.position, ent.cinitial);
    ext = ent.ext;
    saveas (gcf, strcat(DIRECTORY,tagFile,'ConcentracaoInicial',TYPE_TRATAMENT), ext);
end
end

function graphicsCL(ent,fraction, tagTitle, tagFile, ...
    MU,x, A_n_1, Limiar)
% Final concentration with threshold

DIRECTORY = ent.directory; % Output directory
TYPE_TRATAMENT = sprintf('Dot%02ddays%dX%dCi%d', ...
        fraction, ent.days,ent.position, ent.cinitial);
ext = ent.ext;

fig = figure(); set(fig,'visible',ent.view);
plotComLimiar(x,MU*A_n_1,Limiar);
set(gca,'FontSize',12);
xlabel('X (cm)');
ylabel('Concentração (células/cm)');
title(sprintf('Concentração final (%s)',tagTitle));
if ent.save
    saveas (gcf, strcat(DIRECTORY,tagFile,'ConcentracaoLimiar',TYPE_TRATAMENT), ext);
end
end

function graphicsCLComp(ent,fraction, x, AB, clegend)
%Initial Condition Figure
%
fig1 = figure(1); set(fig1,'visible',ent.view);
clf
hold on
for A = AB
    plot (x, A,'LineWidth',1.5);
end
hold off
legend(clegend);
set(gca,'FontSize',12);
ylabel('Concentração (células/cm)');
xlabel('X (cm)');
title('Concentração final');
if ent.save
    DIRECTORY = ent.directory;% Output directory
    TYPE_TRATAMENT = sprintf('Dot%02ddays%dX%dCi%d', ...
        fraction, ent.days,ent.position, ent.cinitial);
    ext = ent.ext;
    saveas (gcf,strcat(DIRECTORY,'ConcentracaoFinal',TYPE_TRATAMENT), ext);
end
end

function graphicsCRM(ent,fraction, tagTitle, tagFile, ...
    MU,x, A_i, A_n_1, R_0, ratio, concentration, Limiar)
% Prints graphs of tumor concentration, radius and maximum concentration.
%
DIRECTORY = ent.directory; % Output directory
TYPE_TRATAMENT = sprintf('Dot%02ddays%dX%dCi%d', ...
        fraction, ent.days,ent.position, ent.cinitial);
ext = ent.ext;
Nt = ent.days * 24;
eixoT = linspace(0, ent.days, Nt);
%initial x final concentration
%clf
fig4 = figure(); set(fig4,'visible',ent.view);
plot (x, A_i,'LineWidth',1.5);
hold on
plot (x, MU*A_n_1,'LineWidth',1.5);
set(gca,'FontSize',12);
xlabel('X (cm)');
ylabel('Concentração (células/cm)');
legend('Concentração inicial', 'Concentração final');
title(sprintf('Condição Inicial %s Final (%s)','{\times}',tagTitle));
hold off
if ent.save
    saveas (gcf, strcat(DIRECTORY,tagFile,'Concentracao',TYPE_TRATAMENT), ext);
end
%Figure evolution of the glioma radius
fig5 = figure(); set(fig5,'visible',ent.view);
plot(eixoT,ratio,'b-','MarkerSize',3, 'LineWidth',2);
hold on
plot(eixoT, R_0*ones(1,Nt),'r');
grid
set(gca,'FontSize',12);
xlabel ('Tempo (days)');
ylabel ('Diâmetro (cm)');
legend('Evolução do Diametro', 'Diametro inicial');
title(sprintf('Diâmetro do Glioma (%s)',tagTitle));
hold off
if ent.save
    saveas (gcf, strcat(DIRECTORY,tagFile,'RaioGlioma',TYPE_TRATAMENT), ext);
end
%maximum concentration
fig6 = figure(); set(fig6,'visible',ent.view);
semilogy(eixoT,MU*concentration,'r-','MarkerSize',3, 'LineWidth',2);
set(gca,'FontSize',12);
xlabel ('tempo (days)');
ylabel ('Concentração máxima (Log(células/cm))');
title(sprintf('Concentração Máxima (%s)',tagTitle));
if ent.save
    saveas (gcf, strcat(DIRECTORY,tagFile,'ConcentracaoMax',TYPE_TRATAMENT), ext);
end
%final concentration with threshold
fig7 = figure(); set(fig7,'visible',ent.view);
plotComLimiar(x,MU*A_n_1,Limiar);
set(gca,'FontSize',12);
xlabel('X (cm)');
ylabel('Concentração (células/cm)');
title(sprintf('Concentração final (%s)',tagTitle));
if ent.save
    saveas (gcf, strcat(DIRECTORY,tagFile,'ConcentracaoLimiar',TYPE_TRATAMENT), ext);
end
end

function graphicsRMComp(ent,fraction, legenda, ...
    MU, R_0, ratioAB, concentrationAB)
% Prints graphs of tumor concentration, radius and maximum concentration.
%
DIRECTORY = ent.directory; % Output directory
TYPE_TRATAMENT = sprintf('Dot%02ddays%dX%dCi%d', ...
        fraction, ent.days,ent.position, ent.cinitial);
Nt = ent.days * 24;
eixoT = linspace(0, ent.days, Nt);
ext = ent.ext;
%Figure evolution of the glioma radius
fig2 = figure(2); set(fig2,'visible',ent.view);
clf
hold on
for ratio = ratioAB
    plot(eixoT,ratio,'MarkerSize',3, 'LineWidth',2);
end
plot(eixoT, R_0*ones(1,Nt),'r');
hold off
grid
set(gca,'FontSize',12);
xlabel ('Tempo (days)');
ylabel ('Diâmetro (cm)');
legend([legenda, {'{R_0}'}],'Location','best');
title('Diâmetro do Glioma');
if ent.save
    saveas (gcf, strcat(DIRECTORY,'compDiametroGlioma',TYPE_TRATAMENT), ext);
end
%maximum concentration
fig3 = figure(3); set(fig3,'visible',ent.view);
for concentration = concentrationAB
    semilogy(eixoT,MU*concentration,'MarkerSize',3, 'LineWidth',2);
    hold on
end
hold off
set(gca,'FontSize',12);
xlabel ('tempo (days)');
ylabel ('Concentração máxima (Log(células/cm))');
legend(legenda);
title('Concentração Máxima');
if ent.save
    saveas (gcf, strcat(DIRECTORY,'compConcentracaoMax',TYPE_TRATAMENT), ext);
end
end
