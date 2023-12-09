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

function [h,a] = plotComLimiar(x,y,limiar)
    lim = y>limiar;
    a=area (x(lim),y(lim),'FaceColor',[0.5,0.5,0.5],'LineStyle',':');
    hold on
    h=plot (x, y,'LineWidth',1.5);
    hold off
    grid
end
