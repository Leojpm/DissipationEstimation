%% Code for the execution of the method in Middleton, Fine, Mackinnon, Alford and Taylor, 2021
% Requires TEOS_10 package


function epsilon = DDCdissipation(S,T,P,lat)

% Create relevant variables
alpha = gsw_alpha(S,T,P);
beta = gsw_beta(S,T,P);
g = gsw_grav(lat,P);
kappa_T = 1.4e-7;
kappa_S = 7.2e-10;
f = gsw_f(lat);

b = g.*alpha.*T-g.*beta.*S; % Create buoyancy field
sp = alpha.*T+beta.*S; % Create spice field
[N2,~] = gsw_Nsquared(S,T,P,lat); % Create buoyancy frequency
N2 = N2(:,1:end-1)+diff(N2,1,2)/2; % Interpolate buoyancy frequency

% Create sorted buoyancy gradient using method of Tseng and Ferziger 2001
[Nn,edges]=histcounts(b,'Normalization','pdf'); 
bgrid = edges(1:end-1);
dzstardb = max(z)*Nn;
dzstardb = interp1(bgrid,dzstardb,b(:));
dzstardb = reshape(dzstardb,size(b));
dzstardb = dzstardb(1:end-1,:)+diff(dzstardb)/2;
dzstardb = dzstardb(:,1:end-1)+diff(dzstardb,1,2)/2;

% Estimate two point correlation R
[X,~] = meshgrid(x,z);
res = 100; % Resolution in buoyancy space
siggrid = min(b(:)):(max(b(:))-min(b(:)))/res:max(b(:));
[Xs,Sig] = meshgrid(x,siggrid);
spgrid = inpaint_nans(griddata(X,inpaint_nans(b),sp,Xs,Sig)); % Interpolate into isopycnal space
Rgrid = ((circshift(spgrid,[0,1])-circshift(spgrid,[0,0])).^2); % Two point correlation
R = griddata(Xs,Sig,Rgrid,X,inpaint_nans(sigma)); % Interpolate back into real space

% Create prefactor in Eq. 10
const1 = -(kappa_T+kappa_S)/2*dzstardb.*sqrt(N.^2./f^2/3+2/3).*modgradb.^2;
const2 = g.*(kappa_T-kappa_S)/2*dzstardb.*sqrt(N.^2./f^2/3+2/3).*modgradb;

% Use numerical method to solve Eq.10 - We iterate, starting from an
% arbitrary initial dissipation, then updating our dissipation estimate
% using Eq. 10 until it converges. Other numerical methods here could be
% used.

epsilon = 1e-5*ones(size(N)); % Initial guess

while nanmean((epsilon(:)-epsilon0(:))./epsilon0(:))<0.01 
    
    koz= 2*pi*N.^(3/2)./epsilon.^(1/2).*f./N; % Create Ozmidov wavenumber in stretched coordinates

    % Create the exponential integral in Eq. 7
    for jj=1:size(rk,1)
        for kk=1:size(rk,2)
            k = logspace(-10,log10(kf(ii,kk)),100);
            int(jj,kk) = trapz(k,k.^(-1).*(1-cos(d(jj,kk)*k))); 
        end
    end
    A = R./int/2; % The prefactor in Eq. 7
    
    epsilon0 = epsilon;
    epsilon = const1 + const2.*sqrt(A/2).*koz; % Eq. 10 estimate for dissipation rate
end


