%% Code for the execution of the method in Middleton, Fine, Mackinnon, Alford and Taylor, 2021
% Requires TEOS_10 package
% Note that this code is for illustrative purposes, so might not apply
% directly to your data

function [epsilon,N] = DDCdissipationV4(S,T,p,lat,lon,kscale)
% The data as input to this code should be gridded vertically in 
% along-track profiles i.e. P should be a single vector for every profile in the matrices T and S.
% The profiles do not need to be a constant distance,
% apart and the method may give more accurate values if you use the true measured
% horizontal distances between profiles as opposed to an artificial 
% horizontal grid. Our final output will be on a half grid, as it uses both
% vertical and horizontal gradients i.e. if your temperature field is an
% (n,m) matrix, your dissipation prediction will be an (n-1,m-1) matrix.

% Inputs
% S : Absolute Salinity on a pressure-time grid, both time and pressure
% spacing may vary.
% T : Conservative Temperature on a pressure-time grid
% p : Pressure (dbar)
% lat : Latitude
% lon : Longitude
% kscale : The exponent for the power spectrum of spice. In Middleton et
% al. (2021) we used k^(-1) (i.e. kscale = -1), following observations from Mackinnon et al.
% (2016) and others. In Fine et al. (2022), we also used k^(-2), as k^(-1)
% overestimated observed dissipation rates.

% Outputs
% epsilon : Turbulent Kinetic Energy Dissipation Rate estimate for averaged
% epsilon between observations based on the Middleton et al. (2021) method.
% N : Buoyancy frequency, N^2 = -dbdz

% Notes:
%  - Dissipation rate has only performed well for low buoyancy Reynolds number
% (Re_b = epsilon/nu/N^2) data. Testing shows values less that 100 work
% well in comparison to measured epsilon: values of 20 have also been used
% as thresholds. I recommend varying this threshold in your data if you are
% making a comparison to microstructure.
%
% - The spacing between observations is important for this method. As shown
% in Middleton et al. (2021), as the spacing increases between
% observations, the skill of the parameterization decreases. Recommended
% spacing is somewhere between 1 km and 10 km. Larger spacing will still
% give an answer, but using a k^(-1) spectrum only will give overinflated
% dissipation rates (using k^(-2) or k^(-3) may give something more
% accurate).
%
% - kscale is limited to values >-3 due to the form of the integral in the
% iterative loop. For values of -3, you should be able to use kscale =
% -2.99.
%
% - The method is best applied to the original sampled data. e.g. if you
% sampled at a variable grid spacing of [1km 1.5km 0.5km 3km] for four
% profiles, then do not interpolate this onto a consistent grid spacing,
% rather use the natural sampling resolution. The same applies to glider
% data. If you can use the spacing between observations derived from the
% original see-saw pattern, you will get better results. For a version of
% this code that works for glider data, feel free to send me an email at
% middleton.leo@gmail.com.
%
% - This method is not a substitute for microstructure observations! It may
% provide insight when no observations are available, but it will not give
% the full picture of turbulent dissipation rate. It can also be used to
% compare to microstructure to understand whether given observations could
% feasibly be explained by double-diffusive processes. However, this
% parameterization does not involve understanding the feedbacks between
% shear-driven turbulence and double-diffusive convection, so it will
% likely underestimate the integrated effects of double diffusion.

%% Create relevant variables
% Create along-track distance x in metres
x(1) = 0;
Ixnan = find(~isnan(lat));
for ii = 2:length(lat)
    deltaLat=lat(ii)*pi/180-lat(max(Ixnan(Ixnan<ii)))*pi/180;
    deltaLon=lon(ii)*pi/180-lon(max(Ixnan(Ixnan<ii)))*pi/180;
    a=sin((deltaLat)/2)^2 + cos(lat(max(Ixnan(Ixnan<ii)))*pi/180)*cos(lat(ii)*pi/180) * sin(deltaLon/2)^2;
    c=2*atan2(sqrt(a),sqrt(1-a));
    dist =6371*c*1000; 
    
    x(ii) = x(max(Ixnan(Ixnan<ii)))+dist;
end

% Make matrix variables
P = repmat(p,length(x),1)'; % Matrix for pressure, make sure the dimensions are the same as T & S
[~,Lat] = meshgrid(p,lat);
[~,Lon] = meshgrid(p,lon);
Z = -gsw_z_from_p(P,Lat');
[~,X] = meshgrid(p,x);
X = X';
Lat = Lat';

% Use the non-linear equation of state to create your buoyancy variable
% i.e. b = g*alpha*T-g*beta*S
alpha = gsw_alpha(S,T,P); % coefficient of thermal expansion
beta = gsw_beta(S,T,P); % coefficient of haline contraction
g = gsw_grav(Lat,P); % Gravity
% b = g.*alpha.*T-g.*beta.*S; % Create buoyancy field
b = -g.*(gsw_rho(S,T,P)-1000)/1000; % using potential density

% Remove overturns by sorting each profile individually
b(isnan(b))=-999;
b = sort(b,1,'descend');
b(b==-999) = nan;

sp = alpha.*T+beta.*S; % Create spice field

dbdx = diff(b')./diff(X'); % horizontal buoyancy gradient
dbdx = dbdx(:,1:end-1)+diff(dbdx,1,2)/2;
dbdz = diff(b)./diff(Z); % vertical buoyancy gradient
dbdz = dbdz';
dbdz = dbdz(1:end-1,:)+diff(dbdz)/2;
modgradb = sqrt(dbdx.^2 + dbdz.^2)'; % Magnitude of the buoyancy gradient

% Create N2

N2 = -dbdz';
N = sqrt(N2);

% Create sorted buoyancy gradient using method of Tseng and Ferziger 2001
[Nn,edges]=histcounts(b,'Normalization','pdf'); % Create pdf of buoyancy
bgrid = edges(1:end-1); % get buoyancy grid
dzstardb = max(Z(:))*Nn; % rescale PDF to get dzstar/db
dzstardb = interp1(bgrid,dzstardb,b(:)); % Interpolate from buoyancy grid to true buoyancy
dzstardb = reshape(dzstardb,size(b)); % Turn into a matrix the same size as b
dzstardb = dzstardb(1:end-1,:)+diff(dzstardb)/2; % Interpolate onto half-grid
dzstardb = dzstardb(:,1:end-1)+diff(dzstardb,1,2)/2;

f = abs(gsw_f(Lat)); % Coriolis parameter

kappa_T = 1.4e-7; % Molecular diffusivity of temperature 
% I recommend using SW_Diffusivity from http://web.mit.edu/seawater/ package
kappa_S = 7.2e-10; % Molecular diffusivity of salt
% I recommend using the equation kappa_S = 1e-9*(0.7+0.036*T) from Vitagliano 1956 

%% Estimate two point correlation R in isopycnal space
% We interpolate each profile into gridded buoyancy space instead of
% gridded depth space. We then calculate a two-point correlation between
% adjacent observations i.e. along-isopycnal, and then we interpolate back
% into gridded depth space.

res = 1000; % Resolution in buoyancy space
bgrid = min(b(:)):(max(b(:))-min(b(:)))/res:max(b(:)); % Create new buoyancy grid, or use the one from above
clearvars spiso R
for ii = 1:length(x) % Interpolate spice onto buoyancy grid
    Inan = find(~isnan(b(:,ii)));
    if isempty(Inan)
        spiso(:,ii) = nan(size(bgrid));
    else
        spiso(:,ii) = interp1(b(Inan,ii),sp(Inan,ii),bgrid);
    end
end
Riso = diff(spiso,1,2).^2; % Calculate the two-point correlation, in horizontal (isopycnal) direction
Riso = Riso(:,1:end-1)+diff(Riso,1,2)/2;
Riso = [nan(size(Riso(:,1))) Riso nan(size(Riso(:,1)))];
db = b(1:end-1,:)+diff(b)/2; % Interpolate buoyancy onto half grid
for ii = 1:length(x)-1 % Interpolate two-point correlation back into real space (now on the half-grid)
    R(ii,:) = interp1(bgrid,Riso(:,ii),db(:,ii));
end
% R = R(:,1:end-1)+diff(R,1,2)/2;

% % Try vertical correlation
% R = diff(sp).^2;
% R = R(:,1:end-1)+diff(R,1,2)/2;
% R = R';

% Put g and f onto half-grid
g = g(1:end-1,:)+diff(g,1,1)/2;
g = g(:,1:end-1)+diff(g,1,2)/2;

f = f(1:end-1,:)+diff(f,1,1)/2;
f = f(:,1:end-1)+diff(f,1,2)/2;

% Create prefactor in Eq. 10
const1 = -(kappa_T+kappa_S)/2*dzstardb.*modgradb.^2; 
const2 = g.*(kappa_T-kappa_S)/2.*dzstardb.*sqrt(N2./f.^2/3+2/3).*modgradb;

% If points are not converging because initially epsilon<0 as const1<0 and |const1|>|const2|
% then try setting const1=0 and adding it back in at the end

% Create horizontal distance field
% D = repmat(diff(x),length(p)-1,1);
D = repmat(diff(p),length(x)-1,1)';

% Use numerical method to solve Eq.10 - We iterate, starting from an
% arbitrary initial dissipation, then updating our dissipation estimate
% using Eq. 10 until it converges. Other numerical methods here could be
% used. 

epsilon0 = 1e-4*ones(size(const2)); % Initial guess 
epsilon = 1e-5*ones(size(const2)); % First guess

while nanmedian(abs((epsilon(:)-epsilon0(:))./epsilon0(:)))>0.1 % We set a threshold for when the prediction at each point has converged
    % Another option is to manually run the loop a few times (e.g. 5) then
    % visualise the whole field epsilon-epsilon0 to check for convergence
    koz= 2*pi*N.^(3/2)./epsilon.^(1/2).*f./N; % Create Ozmidov wavenumber in stretched coordinates

    % Create the exponential integral in Eq. 7
    for jj=1:size(epsilon,1)
        for kk=1:size(epsilon,2)
            k = logspace(-10,log10(koz(jj,kk)),100);
            int(jj,kk) = trapz(k,k.^(kscale).*(1-cos(D(jj,kk)*k))); 
        end
    end
    A = R'./int/2; % The prefactor in Eq. 7
    
    epsilon0 = epsilon;
    epsilon = (const2.*sqrt(A/(3+kscale).*koz.^(3+kscale))); % Eq. 10 estimate for dissipation rate
%     epsilon = (const2.*sqrt(A/2).*koz); % Eq. 10 estimate for dissipation rate


end

