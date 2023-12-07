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
function main_EqRad()
% Example of using the EqRad() function
% The first saves the data in a text file
% The second generates output graphs. 

graphics = {'ci-comp','crm-comp','crm'};
directory = 'img/teste/';
models={'rd', 'rad', 'rd-linear'};
ext_wt = 'txt';
ext_png = 'png';
days = 80;
EqRad(4,'nos',2000,'cinitial',1,'graphics',...
             {'none'}, 'ext', ext_wt, 'view','off','days',days,...
             'models',models,'directory',directory,'wt',true);
EqRad(1,'nos',2000,'cinitial',1,'graphics',...
            graphics, 'ext', ext_png, 'view','off','days',days,...
            'models',models,'directory',directory,'wt',false);
EqRad(4,'nos',2000,'cinitial',1,'graphics',...
            graphics, 'ext', ext_png, 'view','off','days',days,...
            'models',models,'directory',directory,'wt',false);
EqRad(4,'nos',2000,'cinitial',2,'graphics',...
            graphics, 'ext', ext_png, 'view','off','days',days,...
            'models',models,'directory',directory,'wt',false);
EqRad(6,'nos',2000,'cinitial',1,'graphics',graphics,'position',10,...
            'ext', ext_png, 'view','off','days',days,...
            'models',models,'directory',directory,...
            'input', {'D',4e-5,'rho',0.0453,'V',0.0085,'Lx',20,...
            'alpha',0.03,'beta',0.003});
end