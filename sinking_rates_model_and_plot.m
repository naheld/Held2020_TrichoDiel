% Sinking rates for spheres of various density in seawater.
%
% This cannot be the first time such an analysis has been done.
% need to check literature.
%
% Revision History
% 2020-11-03    mvj    Created.
% 2020-11-04    nah    added notes 
% 2020-11-04    nah    added calculation of phytoplankton type cases
% 2020-11-05    nah    make subplot 121 log log, added the main figure
% 2020-12-02    nah    trying to track contour problem in main figure
% 2020-12-04    nah    trying to separate the positive and negatively
% buoyant phytos

clear;

% seawater properties.
NU = 1e-6; % m^2/s
GRAV = 9.81; % m/s^2

% Characteristic dimensions for contours
D = [1e-7 1e-6 1e-5 1e-4 1e-3 1e-2];

% Specific gravities aka density relative to seawater for contours
sigma = [0.2 0.4 0.5 0.8 0.9 0.95 0.98 0.99 1.001 1.002 1.01 1.02 1.05 1.1 1.2 1.5 2 3];

%D = 1e-5;
%sigma = 1.01;

[DD,SS] = meshgrid(D(:),sigma(:));

opt = optimset('MaxFunEvals',1e4,'MaxIter',1e4,'TolX',1e-9);

% read in phytoplankton typecases (density and diameters from the
% literature)
phytoplankton_typecases = readmatrix("phytoplankton_density_typecases_final.csv")
[ppDD, ppSS, ppPP] = meshgrid(phytoplankton_typecases(:,4), phytoplankton_typecases(:,5), phytoplankton_typecases(:,6))
ppsigma = phytoplankton_typecases(:,5)

% above, along with some numerical fits or other means of getting the Cd is enough to solve for the terminal velocity,
% but for general Cd this has to be done numerically.
for mm = 1:size(DD,1)
  for nn = 1:size(DD,2)
    
    %grab the variable of interest
    dd = DD(mm,nn);
    ss = SS(mm,nn);
    
    %starting point J = 1e-4
    [U(mm,nn),J,flag] = fminsearch(@(x) cost(x,dd,ss,NU),1e-4,opt);
    
    %flag is related to reynold's  number - we can fill this in once we
    %have idea of correct Re equations
    if flag ~= 1 % tweak later.
      U(mm,nn) = NaN;
    end

    % store associated Re Fr for contouring later. 
    RE(mm,nn) = Re(U(mm,nn),dd,NU);
    FR(mm,nn) = Fr(U(mm,nn),dd,GRAV);
    
  end
end

% Repeat the calculations for the phytoplankton typecases
for mm = 1:size(ppDD,1)
    
   % one dimensional iterations
   nn = mm
   pp = ppPP(mm,nn);

   %grab the variable of interest
   ppdd = ppDD(mm,nn);
   ppss = ppSS(mm,nn);
    
   %starting point J = 1e-4
   [ppU(pp, mm),J,flag] = fminsearch(@(x) cost(x,ppdd,ppss,NU),1e-4,opt);
    
   %flag is related to reynold's  number - we can fill this in once we
   %have idea of correct Re equations
   if flag ~= 1 % tweak later.
     ppU(pp, mm) = NaN;
   end

   % store associated Re Fr for contouring later. 
   ppRE(pp, mm) = Re(ppU(pp),ppdd,NU);
   ppFR(pp, mm) = Fr(ppU(pp),ppdd,GRAV);
    
  
end

% PLOTTING
% Main text figure

figure(3);
clf reset;

subplot(111); % dimensional results.
[cc,hh] = contour(DD*1e6,SS,U*1e3*60,[-10.^[-6:1:5] 10.^[-6:1:5]],'k--');
clabel(cc,hh);
hold on

% color mapping for the phytoplankton types
cmap = jet(11);
phytoplanktontypes = phytoplankton_typecases(:,6);
[unique_vals, ~, group_idx] = unique(phytoplanktontypes);
col_for_val = cmap(group_idx(:), :); 

scatter(phytoplankton_typecases(:,4)*1e6, phytoplankton_typecases(:,5),120, col_for_val, 'filled');

clabel(cc,hh);
set(gca,'FontSize',14);
xlabel('Diameter (um) - Log scale', 'FontSize', 20);
ylabel({'Specific Gravity (-)'},  'FontSize', 20);
%title('Terminal Velocity (mm/min)',  'FontSize', 20);
set(gca,'xscale','log');
set(gca, 'YDir','reverse');


function J = cost(U,D,sigma,nu)

%(seawater density cancels out, so ignore that in the force balance.
RHOSW = 1000;
 
% drag
re = Re(U,D,nu);
cd = Cd(re);
% extra term -1/2*RHOSW to cancel out seawater density
Fdx = -1/2*RHOSW * U.*abs(U) .* cd .* pi.*(D./2).^2;

% buoyancy
Fbx = 4/3*pi*(D/2).^3 * RHOSW .* (sigma-1); % positive down.

%keyboard

% cost.
J = ((Fbx+Fdx)/1e-24).^2;
%J
%[Fbx Fdx]
%[U re cd Fdx Fbx J]

end


function fr = Fr(U,D,g)
fr = abs(U)./sqrt(g.*D);
end

function re = Re(U,D,nu);
re = abs(U).*D./nu;
end

function cd = Cd(Re)

% Cd calculations
cd = 24/Re ...
    + 2.6*(Re/5.0)/(1+(Re/5.0)^1.52) ...
    + 0.411*(Re/2.63e5)^-7.94/(1+(Re/2.63e5)^-8.00) ...
    + 0.25*(Re/1e6)/(1+(Re/1e6));
end


