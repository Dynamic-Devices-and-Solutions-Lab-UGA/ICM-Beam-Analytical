function [A] = influenceA(a,x,type,BC)
% generate influence matrix A that relates internal moment in a beam to
% external loading as {M} = [A]{F}.
%
% INPUTS:
%
% a - n_loads x 1 - vector that defines load segments. For UDL loading,
%                   a(i) is on the right side of segment i. For LDL loading, a(i) is at the
%                   point of application of p(i).
% x - n_points x 1 - vector defining points along the beam
% type - string - Type of loading ('udl' for uniform distributed loading or
%                 'ldl' for linearly distributed loading)
% BC - string - Type of boundary condition. Options are 'cantilevered'
%               (cantilevered beam), 'simple' (simply supported beam),
%               'clamped' (clamped-clamped beam), and 'clampedsimple' (clamped-simply supported)
%
% OUTPUTS:
%
% A - n_points x n_loads - Influence matrix
%
% EXAMPLES:
%
% A UDL with 4 evenly-spaced control points along a cantilevered, 4m beam.
% Sensors placed at 2.5m and 3.5m
%   [A] = influenceA([1,2,3,4],[2.5,3.5],'udl','cantilevered')
%
% A LDL with 8 evenly-spaced control points along a clamped-clamped, 1m beam.
% Sensors placed at 0.1m, 0.2m, 0.3m, 0.45m, 0.6m, 0.75m, 0.9m, and 0.95m
%   numSens=8;
%   L = 1;
%   a = [0:numSens-1]*L/(numSens-1);
%   sensLocs = [0.1, 0.2, 0.3, 0.45, 0.6, 0.75, 0.9, 0.95];
%   [A] = influenceA(a,sensLocs,'ldl','clamped')

n_loads = length(a);
n_points = length(x);
A = zeros(n_points,n_loads);

switch upper(type)
    case 'UDL'
        if n_loads ~= 1 % If there are more than 1 sections...
            % Define boundary conditions
            Cv = zeros(1,n_loads);
            Cm = Cv;
            switch upper(BC)
                case 'CANTILEVERED'
                    Cv(end) = a(end);
                    Cm(end) = -(a(end)^2)/(2);
                case 'SIMPLE'
                    Cm(1) = (a(1)^2/2);
                    Cm(end) = -(a(end-1)^2/2);
                    for k = 2:length(a)-1
                        Cm(k) = -(-(a(k)^2/2) + (a(k-1)^2/2));
                    end
                    Cv(end) = (1/a(end))*(a(end)^2/2);
                    Cv = Cv - (1/a(end))*(Cm);
                case 'CLAMPED'
                    Cu = Cv; % Initialize displacement and angle BC coefficients
                    Ct = Cv;

                    Cu(1) = (a(1)^4/24);
                    Cu(end) = -(a(end-1)^4/24);
                    for k = 2:length(a)-1
                        Cu(k) = -((a(k-1)^4/24) - (a(k)^4/24));
                    end

                    Ct(1) = -(a(1)^3/6);
                    Ct(end) = (a(end-1)^3/6);
                    for k = 2:length(a)-1
                        Ct(k) = -((a(k)^3/6) - (a(k-1)^3/6));
                    end
                    Cv(end) = (6/a(end)^2)*(a(end)^3/12);
                    Cv = Cv + (6/a(end)^2)*(Ct + (2/a(end))*Cu);

                    Cm(end) = (2/a(end)^2)*(a(end)^4/24);
                    Cm = Cm + (2/a(end)^2)*(-Ct*a(end) - Cu - Cv*(a(end)^3/6));
                case 'CLAMPEDSIMPLE'
                    Cu = Cv; % Initialize displacement and angle BC coefficients
                    Ct = Cv;

                    Cu(1) = (a(1)^4/24);
                    Cu(end) = -(a(end-1)^4/24);
                    for k = 2:length(a)-1
                        Cu(k) = -((a(k-1)^4/24) - (a(k)^4/24));
                    end

                    Ct(1) = -(a(1)^3/6);
                    Ct(end) = (a(end-1)^3/6);
                    for k = 2:length(a)-1
                        Ct(k) = -((a(k)^3/6) - (a(k-1)^3/6));
                    end
                    Cv(end) = ((5/8)*a(end));
                    Cv = Cv + Ct*(3/a(end)^2) + Cu*(3/a(end)^3);

                    Cm(end) = (a(end)^2/2);
                    Cm = Cm - Cv*a(end);
                otherwise
                    error('Invalid support type. Specify a valid type for supportType');
            end

            for i = 1:length(a) % Assemble A by segment
                if i == 1
                    validx = [x >= 0 & x < a(1)];
                else
                    validx = [x >= a(i-1) & x <= a(i)];
                end

                if i ~= length(a) % If x is not in the last segment...
                    A(validx,i) = (a(i)*x(validx) - 0.5*(x(validx).^2 + a(i)^2));
                    A(validx,end) = (((a(end-1)^2)/2)-a(end-1)*x(validx));
                    for k = i+1:length(a)-1
                        A(validx,k) = ((a(k)*x(validx) - ((a(k)^2)/2)) - (a(k-1)*x(validx) - ((a(k-1)^2)/2)));
                    end
                else
                    A(validx,i) = ((-x(validx).^2)/2);
                end
            end
        else % If there's only 1 section...
            % Define boundary conditions
            switch upper(BC)
                case 'CANTILEVERED'
                    Cv = a(1);
                    Cm = -(a(1)^2)/(2);
                case 'SIMPLE'
                    Cm = 0;
                    Cv = a(1)/2;
                case 'CLAMPED'
                    Cv = a(1)/2;
                    Cm = -a(1)^2/12;
                case 'CLAMPEDSIMPLE'
                    Cv = (15*a(1)/24);
                    Cm = -(3*a(1)^2/24);
                otherwise
                    error('Invalid support type. Specify a valid type for supportType');
            end
            A(:,1) = -((x.^2)/2);
        end
    case 'LDL'
        % Define boundary conditions
        Cv = zeros(1,n_loads);
        Cm = Cv;
        
        n = length(a)-1; % number of LDL segments
        mu_first = @(m,x,a) ((-1)/(a(m+1)-a(m)))*(a(m)*(x.^2/2) - (x.^3/6) - a(m)*a(m+1)*x + (a(m+1)^2/2)*x + a(m)*(a(m+1)^2/2) - (a(m+1)^3/3)) + (a(m+1)*x) - (x.^2/2) - (a(m+1)^2/2);
        mu_middle = @(m,x,a) ((-1)/(a(m+1)-a(m)))*(x*((a(m)^2/2) - (a(m)*a(m+1)) + (a(m+1)^2/2)) - (a(m)^3/6) + (a(m)*(a(m+1)^2/2)) - (a(m+1)^3/3)) + (a(m+1)*x) - (x.^2/2) - (a(m+1)^2/2);
        mu_last = @(m,x,a) ((-1)/(a(m+1)-a(m)))*((a(m)^2/2)*x - (a(m)^3/6)) - (x.^2/2);
        % Calculate integrals of mu for boundary condition calculation
        tau_first = @(m,x,a) ((-1)/(a(m+1)-a(m)))*(a(m)*(x.^3/6) - (x.^4/24) - a(m)*a(m+1)*(x.^2/2) + (a(m+1)^2/2)*(x.^2/2) + a(m)*(a(m+1)^2/2)*x - (a(m+1)^3/3).*x) + (a(m+1)*(x.^2/2)) - (x.^3/6) - (a(m+1)^2/2).*x;
        tau_middle = @(m,x,a) ((-1)/(a(m+1)-a(m)))*((x.^2/2)*((a(m)^2/2) - (a(m)*a(m+1)) + (a(m+1)^2/2)) - (a(m)^3/6)*x + (a(m)*(a(m+1)^2/2))*x - (a(m+1)^3/3)*x) + (a(m+1)*(x.^2/2)) - (x.^3/6) - (a(m+1)^2/2)*x;
        tau_last = @(m,x,a) ((-1)/(a(m+1)-a(m)))*((a(m)^2/2)*(x.^2/2) - (a(m)^3/6)*x) - (x.^3/6);
        % ...and 2nd integrals
        nu_first = @(m,x,a) ((-1)/(a(m+1)-a(m)))*(a(m)*(x.^4/24) - (x.^5/120) - a(m)*a(m+1)*(x.^3/6) + (a(m+1)^2/2)*(x.^3/6) + a(m)*(a(m+1)^2/2)*(x.^2/2) - (a(m+1)^3/3).*(x.^2/2)) + (a(m+1)*(x.^3/6)) - (x.^4/24) - (a(m+1)^2/2).*(x.^2/2);
        nu_middle = @(m,x,a) ((-1)/(a(m+1)-a(m)))*((x.^3/6)*((a(m)^2/2) - (a(m)*a(m+1)) + (a(m+1)^2/2)) - (a(m)^3/6)*(x.^2/2) + (a(m)*(a(m+1)^2/2))*(x.^2/2) - (a(m+1)^3/3)*(x.^2/2)) + (a(m+1)*(x.^3/6)) - (x.^4/24) - (a(m+1)^2/2)*(x.^2/2);
        nu_last = @(m,x,a) ((-1)/(a(m+1)-a(m)))*((a(m)^2/2)*(x.^3/6) - (a(m)^3/6)*(x.^2/2)) - (x.^4/24);
        
        switch upper(BC)
            case 'CANTILEVERED'
                if n > 2
                    Cv(n) = a(n+1) - ((-1)/(a(n+1)-a(n)))*(a(n)*a(n+1) - a(n+1)^2/2);
                    Cv(n+1) = -((1)/(a(n+1)-a(n)))*(a(n)*a(n+1) - a(n+1)^2/2);
                    Cm(n) = ((-1)/(a(n+1)-a(n)))*(a(n)*(a(n+1)^2/2) - a(n+1)^3/3) + (-a(n+1)^2/2);
                    Cm(n+1) = ((1)/(a(n+1)-a(n)))*(a(n)*(a(n+1)^2/2) - a(n+1)^3/3);
                elseif n == 2
                    Cv(n) = a(n+1) - ((-1)/(a(n+1)-a(n)))*(a(n)*a(n+1) - a(n+1)^2/2);
                    Cv(n+1) = -((1)/(a(n+1)-a(n)))*(a(n)*a(n+1) - a(n+1)^2/2);
                    Cm(n) = ((-1)/(a(n+1)-a(n)))*(a(n)*(a(n+1)^2/2) - a(n+1)^3/3) + (-a(n+1)^2/2);
                    Cm(n+1) = ((1)/(a(n+1)-a(n)))*(a(n)*(a(n+1)^2/2) - a(n+1)^3/3);
                elseif n ==1
                    Cv(1) = ((-1/(a(2)-a(1)))*(a(2).^2/2 - a(1)*a(2)) + a(2));
                    Cv(2) = ((1/(a(2)-a(1)))*(a(2).^2/2 - a(1)*a(2)));
                    Cm(1) = ((-1/(a(2)-a(1)))*(a(2).^3/6 - a(1)*a(2).^2/2) + a(2)^2/2);
                    Cm(2) = ((1/(a(2)-a(1)))*(a(2).^3/6 - a(1)*a(2).^2/2));
                    Cm = Cm - Cv*a(2);
                end
            case 'SIMPLE'
                if n > 2
                    Cm(1) = -mu_first(1,0,a);
                    Cm(2) = -(mu_middle(2,0,a) - mu_first(1,0,a));
                    Cm(n) = -(mu_last(n,0,a) - mu_middle(n-1,0,a));
                    Cm(n+1) = mu_last(n,0,a);
                    for k = 3:n-1
                        Cm(k) = -(mu_middle(k,0,a) - mu_middle(k-1,0,a));
                    end

                    Cv(n) = (1/(a(n+1)-a(n)))*(-a(n+1)^3/6 + a(n)*(a(n+1)^2/2)) + (a(n+1)^2/2);
                    Cv(n+1) = (-1/(a(n+1)-a(n)))*(-a(n+1)^3/6 + a(n)*(a(n+1)^2/2));
                    Cv = (1/a(n+1))*(Cv-Cm);
                elseif n == 2
                    Cm(1) = -mu_first(1,0,a);
                    Cm(2) = -(mu_last(2,0,a)-mu_first(1,0,a));
                    Cm(3) = mu_last(2,0,a);
                    Cv(1) = 0;
                    Cv(2) = (1/a(3))*(((1/(a(3)-a(2)))*(-a(3)^3/6 + a(2)*a(3)^2/2) + a(3)^2/2));
                    Cv(3) = (1/a(3))*(((-1/(a(3)-a(2)))*(-a(3)^3/6 + a(2)*a(3)^2/2)));
                    Cv = Cv + (1/a(3))*(-Cm);
                elseif n ==1
                    Cm = 0;
                    Cv(1) = (1/a(2))*(((-1/(a(2)-a(1)))*(a(2).^3/6 - a(1)*a(2).^2/2) + a(2)^2/2));
                    Cv(2) = (1/a(2))*(((1/(a(2)-a(1)))*(a(2).^3/6 - a(1)*a(2).^2/2)));
                end
            case 'CLAMPED'
                Cu = Cv;
                Ct = Cv;
                if n > 2
                    Cu(1) = (nu_first(1,a(2),a) - tau_first(1,a(2),a)*a(2));
                    Cu(2) = (-nu_first(2,a(2),a) + nu_middle(2,a(2),a) - nu_first(1,a(2),a) + tau_first(2,a(2),a)*a(2) + tau_first(1,a(2),a)*a(2) - tau_middle(2,a(2),a)*a(2) + nu_first(2,a(3),a) - tau_first(2,a(3),a)*a(3));
                    Cu(n) = (nu_first(n-1,a(n-1),a) - nu_middle(n-1,a(n-1),a) + tau_middle(n-1,a(n-1),a)*a(n-1) - tau_first(n-1,a(n-1),a)*a(n-1) + (1/(a(n+1)-a(n)))*(a(n)^5/30) + (a(n)^4/24) + nu_last(n,a(n),a) - nu_first(n-1,a(n),a) + (-1/(a(n+1)-a(n)))*(a(n)^5/8) - (a(n)^4/6) - tau_last(n,a(n),a)*a(n) + tau_first(n-1,a(n),a)*a(n));
                    Cu(n+1) = ((-1/(a(n+1)-a(n)))*(a(n)^5/30) - nu_last(n,a(n),a) - (a(n)^4/24) + (1/(a(n+1)-a(n)))*(a(n)^5/8) + tau_last(n,a(n),a)*a(n) + (a(n)^4/6));
                         for k = 3:n-1
                             Cu(k) = (nu_first(k-1,a(k-1),a) - nu_middle(k-1,a(k-1),a) - tau_first(k-1,a(k-1),a)*a(k-1) + tau_middle(k-1,a(k-1),a)*a(k-1) - nu_first(k-1,a(k),a) + tau_first(k-1,a(k),a)*a(k) - nu_first(k,a(k),a) + nu_middle(k,a(k),a) + tau_first(k,a(k),a)*a(k) - tau_middle(k,a(k),a)*a(k) + nu_first(k,a(k+1),a) - tau_first(k,a(k+1),a)*a(k+1));
                         end
                    Ct(1) = (tau_first(1,a(2),a));
                    Ct(2) = (-tau_first(2,a(2),a) - tau_first(1,a(2),a) + tau_middle(2,a(2),a) + tau_first(2,a(3),a));
                    Ct(n) = (-tau_middle(n-1,a(n-1),a) + tau_first(n-1,a(n-1),a) + (1/(a(n+1)-a(n)))*(a(n)^4/8) + (a(n)^3/6) + tau_last(n,a(n),a) - tau_first(n-1,a(n),a));
                    Ct(n+1) = ((-1/(a(n+1)-a(n)))*(a(n)^4/8) - tau_last(n,a(n),a) - (a(n)^3/6));
                         for k = 3:n-1
                             Ct(k) = (tau_first(k-1,a(k-1),a) - tau_middle(k-1,a(k-1),a) - tau_first(k-1,a(k),a) - tau_first(k,a(k),a) + tau_middle(k,a(k),a) + tau_first(k,a(k+1),a));
                         end

                    Cm(n) = (6/(a(n+1)^2))*(((1/(a(n+1)-a(n)))*((-a(n+1)^5/120) + a(n)*(a(n+1)^4/24)) + (a(n+1)^4/24))) - (2/(a(n+1)))*(((1/(a(n+1)-a(n)))*((-a(n+1)^4/24) + a(n)*(a(n+1)^3/6)) + (a(n+1)^3/6)));
                    Cm(n+1) = (6/(a(n+1)^2))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^5/120) + a(n)*(a(n+1)^4/24)))) - (2/(a(n+1)))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^4/24) + a(n)*(a(n+1)^3/6))));
                    Cm = Cm + (6/(a(n+1)^2))*(-Ct*a(n+1) - Cu) - (2/(a(n+1)))*(-Ct);
                    Cv(n) = (2/(a(n+1)^2))*(((1/(a(n+1)-a(n)))*((-a(n+1)^4/24) + a(n)*(a(n+1)^3/6)) + (a(n+1)^3/6)));
                    Cv(n+1) = (2/(a(n+1)^2))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^4/24) + a(n)*(a(n+1)^3/6))));
                    Cv = Cv + (2/(a(n+1)^2))*(-Cm*a(n+1) - Ct);
                elseif n == 2
                    Cu(1) = -(-nu_first(1,a(2),a) + tau_first(1,a(2),a)*a(2));
                    Cu(2) = -((-1/(a(3)-a(2)))*(a(2)^5/30) - a(2)^4/24 - nu_last(2,a(2),a) + nu_first(1,a(2),a) + (1/(a(3)-a(2)))*(a(2)^5/8) + a(2)^4/6 + tau_last(2,a(2),a)*a(2) - tau_first(1,a(2),a)*a(2));
                    Cu(3) = -((1/(a(3)-a(2)))*(a(2)^5/30) + nu_last(2,a(2),a) + a(2)^4/24 + (-1/(a(3)-a(2)))*(a(2)^5/8) - tau_last(2,a(2),a)*a(2) - a(2)^4/6);
                    
                    Ct(1) = -(-tau_first(1,a(2),a));
                    Ct(2) = -((-1/(a(3)-a(2)))*(a(2)^4/8) - a(2)^3/6 - tau_last(2,a(2),a) + tau_first(1,a(2),a));
                    Ct(3) = -((1/(a(3)-a(2)))*(a(2)^4/8) + tau_last(2,a(2),a) + a(2)^3/6);
                    
                    Cm(n) = (6/(a(n+1)^2))*(((1/(a(n+1)-a(n)))*((-a(n+1)^5/120) + a(n)*(a(n+1)^4/24)) + (a(n+1)^4/24))) - (2/(a(n+1)))*(((1/(a(n+1)-a(n)))*((-a(n+1)^4/24) + a(n)*(a(n+1)^3/6)) + (a(n+1)^3/6)));
                    Cm(n+1) = (6/(a(n+1)^2))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^5/120) + a(n)*(a(n+1)^4/24)))) - (2/(a(n+1)))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^4/24) + a(n)*(a(n+1)^3/6))));
                    Cm = Cm + (6/(a(n+1)^2))*(-Ct*a(n+1) - Cu) - (2/(a(n+1)))*(-Ct);
                    
                    Cv(n) = (2/(a(n+1)^2))*(((1/(a(n+1)-a(n)))*((-a(n+1)^4/24) + a(n)*(a(n+1)^3/6)) + (a(n+1)^3/6)));
                    Cv(n+1) = (2/(a(n+1)^2))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^4/24) + a(n)*(a(n+1)^3/6))));
                    Cv = Cv + (2/(a(n+1)^2))*(-Cm*a(n+1) - Ct);
                elseif n ==1
                    Cv(1) = (-6/a(2)^2)*(((1/(a(2)-a(1)))*(a(2)^4/24 - a(1)*a(2)^3/6) - a(2)^3/6) + (-2/a(2))*(((1/(a(2)-a(1)))*(a(2)^5/120 - a(1)*a(2)^4/24) - a(2)^4/24)));
                    Cv(2) = (-6/a(2)^2)*(((-1/(a(2)-a(1)))*(a(2)^4/24 - a(1)*a(2)^3/6)) + (-2/a(2))*(((-1/(a(2)-a(1)))*(a(2)^5/120 - a(1)*a(2)^4/24))));
                    Cm(1) = (-2/a(2)^2)*(((1/(a(2)-a(1)))*(a(2)^5/120 - a(1)*a(2)^4/24) - a(2)^4/24));
                    Cm(2) = (-2/a(2)^2)*(((-1/(a(2)-a(1)))*(a(2)^5/120 - a(1)*a(2)^4/24)));
                    Cm = Cm + (-2/a(2)^2)*(Cv*(a(2)^3/6));
                    
                end
                
            case 'CLAMPEDSIMPLE'
                Cu = Cv;
                Ct = Cv;
                if n > 2
                    Cu(1) = (nu_first(1,a(2),a) - tau_first(1,a(2),a)*a(2));
                    Cu(2) = (-nu_first(2,a(2),a) + nu_middle(2,a(2),a) - nu_first(1,a(2),a) + tau_first(2,a(2),a)*a(2) + tau_first(1,a(2),a)*a(2) - tau_middle(2,a(2),a)*a(2) + nu_first(2,a(3),a) - tau_first(2,a(3),a)*a(3));
                    Cu(n) = (nu_first(n-1,a(n-1),a) - nu_middle(n-1,a(n-1),a) + tau_middle(n-1,a(n-1),a)*a(n-1) - tau_first(n-1,a(n-1),a)*a(n-1) + (1/(a(n+1)-a(n)))*(a(n)^5/30) + (a(n)^4/24) + nu_last(n,a(n),a) - nu_first(n-1,a(n),a) + (-1/(a(n+1)-a(n)))*(a(n)^5/8) - (a(n)^4/6) - tau_last(n,a(n),a)*a(n) + tau_first(n-1,a(n),a)*a(n));
                    Cu(n+1) = ((-1/(a(n+1)-a(n)))*(a(n)^5/30) - nu_last(n,a(n),a) - (a(n)^4/24) + (1/(a(n+1)-a(n)))*(a(n)^5/8) + tau_last(n,a(n),a)*a(n) + (a(n)^4/6));
                         for k = 3:n-1
                             Cu(k) = (nu_first(k-1,a(k-1),a) - nu_middle(k-1,a(k-1),a) - tau_first(k-1,a(k-1),a)*a(k-1) + tau_middle(k-1,a(k-1),a)*a(k-1) - nu_first(k-1,a(k),a) + tau_first(k-1,a(k),a)*a(k) - nu_first(k,a(k),a) + nu_middle(k,a(k),a) + tau_first(k,a(k),a)*a(k) - tau_middle(k,a(k),a)*a(k) + nu_first(k,a(k+1),a) - tau_first(k,a(k+1),a)*a(k+1));
                         end
                    Ct(1) = (tau_first(1,a(2),a));
                    Ct(2) = (-tau_first(2,a(2),a) - tau_first(1,a(2),a) + tau_middle(2,a(2),a) + tau_first(2,a(3),a));
                    Ct(n) = (-tau_middle(n-1,a(n-1),a) + tau_first(n-1,a(n-1),a) + (1/(a(n+1)-a(n)))*(a(n)^4/8) + (a(n)^3/6) + tau_last(n,a(n),a) - tau_first(n-1,a(n),a));
                    Ct(n+1) = ((-1/(a(n+1)-a(n)))*(a(n)^4/8) - tau_last(n,a(n),a) - (a(n)^3/6));
                         for k = 3:n-1
                             Ct(k) = (tau_first(k-1,a(k-1),a) - tau_middle(k-1,a(k-1),a) - tau_first(k-1,a(k),a) - tau_first(k,a(k),a) + tau_middle(k,a(k),a) + tau_first(k,a(k+1),a));
                         end

                    Cm(n) = (3/(a(n+1)^2))*(((1/(a(n+1)-a(n)))*((-a(n+1)^5/120) + a(n)*(a(n+1)^4/24)) + (a(n+1)^4/24)) - (a(n+1)^2/6)*(((1/(a(n+1)-a(n)))*((-a(n+1)^3/6) + a(n)*(a(n+1)^2/2)) + (a(n+1)^2/2))));
                    Cm(n+1) = (3/(a(n+1)^2))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^5/120) + a(n)*(a(n+1)^4/24))) - (a(n+1)^2/6)*(((-1/(a(n+1)-a(n)))*((-a(n+1)^3/6) + a(n)*(a(n+1)^2/2)))));
                    Cm = Cm + (3/(a(n+1)^2))*(-Ct*a(n+1) - Cu);

                    Cv(n) = (1/a(n+1))*(((1/(a(n+1)-a(n)))*((-a(n+1)^3/6) + a(n)*(a(n+1)^2/2)) + (a(n+1)^2/2)));
                    Cv(n+1) = (1/a(n+1))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^3/6) + a(n)*(a(n+1)^2/2))));
                    Cv = Cv + (1/a(n+1))*(-Cm);
                elseif n == 2
                    Cu(1) = -(-nu_first(1,a(2),a) + tau_first(1,a(2),a)*a(2));
                    Cu(2) = -((-1/(a(3)-a(2)))*(a(2)^5/30) - a(2)^4/24 - nu_last(2,a(2),a) + nu_first(1,a(2),a) + (1/(a(3)-a(2)))*(a(2)^5/8) + a(2)^4/6 + tau_last(2,a(2),a)*a(2) - tau_first(1,a(2),a)*a(2));
                    Cu(3) = -((1/(a(3)-a(2)))*(a(2)^5/30) + nu_last(2,a(2),a) + a(2)^4/24 + (-1/(a(3)-a(2)))*(a(2)^5/8) - tau_last(2,a(2),a)*a(2) - a(2)^4/6);
                    
                    Ct(1) = -(-tau_first(1,a(2),a));
                    Ct(2) = -((-1/(a(3)-a(2)))*(a(2)^4/8) - a(2)^3/6 - tau_last(2,a(2),a) + tau_first(1,a(2),a));
                    Ct(3) = -((1/(a(3)-a(2)))*(a(2)^4/8) + tau_last(2,a(2),a) + a(2)^3/6);
                    
                    Cm(n) = (3/(a(n+1)^2))*(((1/(a(n+1)-a(n)))*((-a(n+1)^5/120) + a(n)*(a(n+1)^4/24)) + (a(n+1)^4/24)) - (a(n+1)^2/6)*(((1/(a(n+1)-a(n)))*((-a(n+1)^3/6) + a(n)*(a(n+1)^2/2)) + (a(n+1)^2/2))));
                    Cm(n+1) = (3/(a(n+1)^2))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^5/120) + a(n)*(a(n+1)^4/24))) - (a(n+1)^2/6)*(((-1/(a(n+1)-a(n)))*((-a(n+1)^3/6) + a(n)*(a(n+1)^2/2)))));
                    Cm = Cm + (3/(a(n+1)^2))*(-Ct*a(n+1) - Cu);

                    Cv(n) = (1/a(n+1))*(((1/(a(n+1)-a(n)))*((-a(n+1)^3/6) + a(n)*(a(n+1)^2/2)) + (a(n+1)^2/2)));
                    Cv(n+1) = (1/a(n+1))*(((-1/(a(n+1)-a(n)))*((-a(n+1)^3/6) + a(n)*(a(n+1)^2/2))));
                    Cv = Cv + (1/a(n+1))*(-Cm);
                elseif n ==1
                    Cv(1) = (-3/(2*a(2)))*(((1/(a(2)-a(1)))*(a(2)^3/6 - a(1)*a(2)^2/2) - a(2)^2/2) + (-2/a(2)^2)*(((1/(a(2)-a(1)))*(a(2)^5/120 - a(1)*a(2)^4/24) - a(2)^4/24)));
                    Cv(2) = (-3/(2*a(2)))*(((-1/(a(2)-a(1)))*(a(2)^3/6 - a(1)*a(2)^2/2)) + (-2/a(2)^2)*(((-1/(a(2)-a(1)))*(a(2)^5/120 - a(1)*a(2)^4/24))));
                    Cm(1) = (-2/a(2)^2)*(((1/(a(2)-a(1)))*(a(2)^5/120 - a(1)*a(2)^4/24) - a(2)^4/24));
                    Cm(2) = (-2/a(2)^2)*(((-1/(a(2)-a(1)))*(a(2)^5/120 - a(1)*a(2)^4/24)));
                    Cm = Cm + (-2/a(2)^2)*(Cv*(a(2)^3/6));
                end
            otherwise
                error('Invalid support type. Specify a valid type for supportType');
        end
        
        for i = 1:n
            validx = [x >= a(i) & x <= a(i+1)];

            if i < n-1 % If x is not in the last segment...
                A(validx,i) = mu_first(i,x(validx),a);
                A(validx,i+1) = (mu_middle(i+1,x(validx),a) - mu_first(i,x(validx),a));
                A(validx,n) = (mu_last(n,x(validx),a) - mu_middle(n-1,x(validx),a));
                A(validx,n+1) = -(mu_last(n,x(validx),a) + (x(validx).^2/2));
                for k = i+2:n-1
                    A(validx,k) = (mu_middle(k,x(validx),a) - mu_middle(k-1,x(validx),a));
                end
            elseif i == n-1
                A(validx,i) = mu_first(i,x(validx),a);
                A(validx,i+1) = (mu_last(i+1,x(validx),a) - mu_first(i,x(validx),a));
                A(validx,i+2) = -(mu_last(n,x(validx),a) + (x(validx).^2/2));
            elseif i == n
                A(validx,i) = ((-1/(a(i+1)-a(i)))*(a(i)*(x(validx).^2/2) - x(validx).^3/6)-x(validx).^2/2);
                A(validx,i+1) = ((1/(a(i+1)-a(i)))*(a(i)*(x(validx).^2/2) - x(validx).^3/6));
            end
        end
    otherwise
        error('Invalid Load type. Specify a valid type for loadType');
end
for i = 1:n_points % Add boundary conditions
    A(i,:) = A(i,:) + Cv*x(i) + Cm;
end
end