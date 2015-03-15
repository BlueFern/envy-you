classdef WallMechanics < handle
    % The 'WallMechanics' code contains the following sections of the model:
    %   The Contraction and Mechanical Models.
    %   Please refer to the relevient sections in the documentation for 
    %   full information on the equations and variable names.
    properties
        params
        u0
        index
        n_out
        idx_out
        enabled
    end
    methods
        function self = WallMechanics(varargin)
            self.params = parse_inputs(varargin{:});
            self.index = indices();
            self.u0 = initial_conditions(self.index);
            self.enabled = true(size(self.u0));
            [self.idx_out, self.n_out] = output_indices();
        end
        function [du, varargout] = rhs(self, t, u, Ca_i)
            % Initalise inputs and parameters
            idx = self.index;
            p = self.params;
            % Make the output function compute these so as to avoid
            % duplicating equations
            [R, h] = self.shared(t, u);
            
            Mp = u(idx.Mp, :);
            AMp = u(idx.AMp, :);
            AM = u(idx.AM, :);
            
            %% Contraction Equations
            K_1 = p.gamma_cross * Ca_i.^3;
            K_6 = K_1;
            M = 1 - AM - AMp - Mp;
            du(idx.Mp, :) = p.K_4 * AMp + K_1 .* M - (p.K_2 + p.K_3) * Mp;
            du(idx.AMp, :) = p.K_3 * Mp + K_6 .* AM - (p.K_4 + p.K_5) * AMp;
            du(idx.AM, :) = p.K_5 * AMp - (p.K_7 + K_6) .* AM;
            
            % Mechanical Equations
            F_r = AMp + AM;
            E = p.E_passive + F_r * (p.E_active - p.E_passive);
            R_0 = p.R_0_passive + F_r * (p.alpha - 1) * p.R_0_passive;
            du(idx.R, :) = p.R_0_passive / p.eta * ( R * p.P_T ./ h - ...
                E .* (R - R_0) ./ R_0);
            
            du = bsxfun(@times, self.enabled, du);
            
            if nargout == 2
               Uout = zeros(self.n_out, size(u, 2));
               Uout(self.idx_out.M, :) = M;
               Uout(self.idx_out.F_r, :) = F_r;
               varargout{1} = Uout;
            end
        end
        function [R, h] = shared(self, ~, u)
            R = u(self.index.R, :);
            h = 0.1 * R;
        end     
        function names = varnames(self)
            names = [fieldnames(self.index); fieldnames(self.idx_out)];
        end
    end 
end

function idx = indices()
% Index of parameters needing inital conditions 
idx.Mp = 1;
idx.AMp = 2;
idx.AM = 3;
idx.R = 4;
end

function [idx, n] = output_indices()
% Index of all other output parameters
idx.M = 1;
idx.F_r = 2;
n = numel(fieldnames(idx));
end

function params = parse_inputs(varargin)
parser = inputParser();
% Contraction Equation Constants
parser.addParameter('K_2', 0.5); % s^-1
parser.addParameter('K_3', 0.4); % s^-1
parser.addParameter('K_4', 0.1); % s^-1
parser.addParameter('K_5', 0.5); % s^-1
parser.addParameter('K_7', 0.1); % s^-1
parser.addParameter('gamma_cross', 17); %uM^-3 s^-1
% Mechanical Equation Constants
parser.addParameter('eta', 1e4); %Pa s
parser.addParameter('R_0_passive', 20e-6); % m
parser.addParameter('h_0_passive', 3e-6); % m
parser.addParameter('P_T', 4000); % Pa
parser.addParameter('E_passive', 66e3); % Pa
parser.addParameter('E_active', 233e3); % Pa
parser.addParameter('alpha', 0.6); % [-]

parser.parse(varargin{:});
params = parser.Results;

end
function u0 = initial_conditions(idx)
u0 = zeros(length(fieldnames(idx)), 1);
% Inital estimations of parameters from experimental data
u0(idx.Mp) = .25;
u0(idx.AMp) = .25;
u0(idx.AM) = .25;
u0(idx.R) = 15e-6;
end
