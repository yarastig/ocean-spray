function spray24ghb
%Solving the full thermodynamic equations including the spray source terms
%at the wavecrest level. Boundary conditions are spesified at z1<zw.
%Ventilation effect is included.
%Authors: Sergey A. Suslov & Yevgenii Rastigejev
%Last modified: 04.07.26
%==========================================================================
close all;
disp('-------------------------------------------------------------------')
global N np alp alpe alpeps alpv alpT sigma s0 r0 a0 D0 sig0t Cb1 Cb2
global pi0 pi1 pi2 pi3 pi5 pi6 pi7 pi8 pi9 pi12 pi13 pi17 l U tau cpd cpv
global nu g us kd q0 Ceps1 Ceps2 Ceps3 kp
global dr ri pii s0i a0i
global z0 z1 zq zinf u1log qs thetas T0 P0 thetainf qinf
global Sols1 Sols2 Sols3 Sols4 Sols5 Sols6 Sols7 Sols8 Sols9 Sols10
global Sols11 Sols12 Sols13 Sols14 Sols15 Sols16 Sols17
global kas C0is C1is ais dis kis xiris
global options odeoptions tola tolit nitmax zint zplt is0 ithq


%**** Reference case (no spray):
ithq=0;% ithq=0: thetainf and qinf are computed for the first value of s00
       % ithq=1: thetainf and qinf are read from the initial guess file
%**** Roughness length:
iz0=0; % iz0=0 Charnock roughness length is used
       % iz0=1 Large & Pond (1981) roughness length is used
%**** Initialisation:
igs=1; % igs=0 initial guess is set by this code
       % igs=1 initial guess is read from file fin
fin='Ns1.0e-09us4.0R060.dat';
%==========================================================================
s00=10^(-9);
ns0=length(s00);% array of bin concentrations at the wavecrest level
us=4.;         % m/s, friction velocity
qs=0.0004;      % parameter entering the boundary condition for q
thetas=0.2;     % parameter entering the boundary condition for Ta
U=0/us;         % scaled spray injection speed
%==========================================================================
%%%%% Fixed physical parameters:
zwt=5;           % m, wavecrest level
h=1;             % m, the half thickness of the spray production layer
g=9.8;           % m/s^2, gravity acceleration
nu=1.568*10^(-5);% m^2/s, kinematic viscosity of dry air at T=275K
rhow=1020;       % kg/m^3, density of seawater
Rd=287;          % J/kg/K, gas constant for dry air
Rv=461.51;       % J/kg/K, gas constant for water vapour
cpv=1870;        % J/kg/K, heat capacity of water vapour
cw=4218;         % J/kg/K, heat capacity of water
cpd=1005;        % J/kg/K, heat capacity of dry air
l0=2.45*10^6;    % J/kg, latent heat for T=300K
kd=2/7;          % parameter in the definition of theta
l=l0/l0;         % non-dimensional latent heat (assumed to be constant)

%**** Turbulence and model fitting parameters:
Ceps1=1.44;      % k-epsilon model parameter
Ceps2=1.92;      % k-epsilon model parameter
Ceps3=1;         % k-epsilon model parameter
alp=0.3;         % k-epsilon model parameter
alpe=1;          % k-epsilon model parameter
alpeps=0.77;     % k-epsilon model parameter
alpv=1;          % vapour diffusion constant
alpT=alpv;       % thermal diffusion constant
sig0t=0.67;      % constant for relaxation model
Cb1=0.45;        % constant entering definition of C0
Cb2=1.35;        % constant entering definition of C0
kp=sqrt((Ceps2-Ceps1)*alp/alpeps); % von Karman constant
tau=sqrt(2*alp/alpe)/kp;           % constant entering bc for ea
z0=roughl(us,iz0,g,kp)/zwt; % the nondimensional lower edge of the
                            % turbulent boundary layer
uwlog=log(1/z0)/kp;         % nondimensional logarithmic velocity at the
                            % wavecrest level
Rr=us*z0*zwt/nu;               % roughness Reynolds number
zq=1.e-4*min(1.1,0.55*Rr^(-0.6))/zwt; % the nondimensional thermal
                                      % roughness length
%**** Reference values:
T0=300;          % K, reference temperature
P0=101300;       % Pa, reference pressure of dry air
rhod0=P0/Rd/T0;  % kg/m^3, reference density of dry air
D0=-1.92e1*(1-7.00e-3*T0+1.07e-5*T0^2)/P0;  % m^2/s, reference diffusion
                                            % coefficient of vapour in air
K0=4.64e-2*(1-4.69e-3*T0+1.08e-5*T0^2);     % J/m/s/K, reference thermal
                                            % conductivity of vapour
pvs0=611*exp(53.49-6808/T0-5.09*log(T0));   % Pa, reference vapour pressure
q0=Rd/Rv/P0*pvs0;                           % reference humidity 
sigma=rhow/rhod0; % water/air density ratio
%**** Spray independent nondimensional parameters:
pi0=zwt/2/h         
pi8=g*zwt/Rd/T0
pi9=(Rv/Rd-1)*q0
pi13=alpT*g*zwt/us^2
%--------------------------------------------------------------------------
%**** Computational parameters:
np=5001;     % the number of discretization points in the layer
zinf=200;    % the top edge of the nondimensional computational domain
z1=4*z0;     % the lower edge of the nondimensional computational domain,
             % 0<z1<1
u1log=log(z1/z0)/kp;        % nondimensional logarithmic velocity at the
                            % lower edge of the computational domain
zint=linspace(z1,zinf,np)'; % computational grid
N=1;         % the number of droplet size bins
Nint=1000;   % the number of points for distribution interpolation
tola=1.e-10; % absolute tolerance for terminal speed calculations
options=optimset('TolFun',tola,'Display','off');
tolr1=1.e-3; % relative tolerance in ODE solver
tola1=1.e-4; % absolute tolerance in ODE solver
nmx=75000;   % the maximum number of points
odeoptions=bvpset('RelTol',tolr1,'AbsTol',tola1,'NMax',nmx);
tolit=.05;   %iteration error
nitmax=20;   %the maximum number of iterations
%**** Plotting parameters:
fsz=16; lw=2;
zplt=20;     % z plotting range
%***** Spectral distributions:
Rsmin=0.000010;   % m, left cutoff of the known spectrum 
N=1;
tp='N'; % output data file prefix for single bin spray distribution
if iz0==1
    tp='NLP';
end
Rmin=0.000060; % m, droplet radius 
Rmax=0.000060; % m, droplet radius 
rit=Rmin       % m, droplet radius
drt=Rmin/2;    % m, half-width of the bin
s0i=1;         % spray concentartion
f=@(grb)termv(grb,sigma,nu,Rmin,g);
a=fsolve(f,1,options);
%**** Equivalent monodisperse spray terminal velocity:
a0=a; r0=Rmin;
a0i=a; pii=1;
ri=rit/r0;                      % non-dimensional droplet radii
dr=drt/r0;                      % non-dimensional bin half-width
%**** Spray-dependent nondimensional parameters:      
pi1=a0/us      
pi2=3*D0*q0*zwt/sigma/us/r0^2
pi6=3*K0*zwt/us/r0^2/cw/rhow
pi7=3*D0*q0*zwt*l0/cw/sigma/us/r0^2/T0
pi17=a0*us/g/zwt

%**************************** Initialisation ******************************
ig=zeros(np,11+6*N); xirig=zeros(np,N);
xiri(1:np,1:N)=.1; % initial guess for xiri
if igs==0
    for i=1:np
        ig(i,:)=iguess(zint(i))';
    end
else
    fid=fopen(fin);
        fgetl(fid);
        arr=fscanf(fid,'%e',6);
        Nig=fix(arr(1)), npig=fix(arr(2))
        if ithq==1
            thetainf=arr(5), qinf=arr(6)
        end
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid);
        fscanf(fid,'%e',6);
        fgetl(fid); fgetl(fid); 
        zig=zeros(npig,1); igig=zeros(npig,11+6*Nig);
        for i=1:npig
            arr=fscanf(fid,'%e',7);
            zig(i)=arr(1); igig(i,1:6)=arr(2:7);
        end
        for i=1:6
        	ig(:,i)=spline(zig,igig(:,i),zint);
        end
        fgetl(fid); fgetl(fid); 
        for i=1:npig
            arr=fscanf(fid,'%e',7);
            igig(i,7:11)=arr(2:6);
        end
        for i=1:5
        	ig(:,6+i)=spline(zig,igig(:,6+i),zint);
        end
        xirigig=zeros(npig,Nig);
        for j=1:Nig
            fgetl(fid); fgetl(fid);
            fgetl(fid);
            fscanf(fid,'%e',3);
            fgetl(fid); fgetl(fid);
            for i=1:npig
                arr=fscanf(fid,'%e',7);
                igig(i,0*Nig+11+j)=arr(2);
                igig(i,1*Nig+11+j)=arr(3);
                igig(i,2*Nig+11+j)=arr(4);
                igig(i,3*Nig+11+j)=arr(5);
                igig(i,4*Nig+11+j)=arr(6);
                igig(i,5*Nig+11+j)=arr(7);
            end
            fgetl(fid); fgetl(fid);
            for i=1:npig
                arr=fscanf(fid,'%e',7);
                xirigig(i,j)=arr(7);
            end
        end
    fclose(fid);
    mnn=min(N,Nig);
    for j=1:mnn
        ig(:,0*N+11+j)=spline(zig,igig(:,0*Nig+11+j),zint);
        ig(:,1*N+11+j)=spline(zig,igig(:,1*Nig+11+j),zint);
        ig(:,2*N+11+j)=spline(zig,igig(:,2*Nig+11+j),zint);
        ig(:,3*N+11+j)=spline(zig,igig(:,3*Nig+11+j),zint);
        ig(:,4*N+11+j)=spline(zig,igig(:,4*Nig+11+j),zint);
        ig(:,5*N+11+j)=spline(zig,igig(:,5*Nig+11+j),zint);
        xirig(:,j)=spline(zig,xirigig(:,j),zint);
    end
end
Sols1=spline(zint,ig(:,1));
Sols2=spline(zint,ig(:,2));
Sols3=spline(zint,ig(:,3));
Sols4=spline(zint,ig(:,4));
Sols5=spline(zint,ig(:,5));
Sols6=spline(zint,ig(:,6));
Sols7=spline(zint,ig(:,7));
Sols8=spline(zint,ig(:,8));
Sols9=spline(zint,ig(:,9));
Sols10=spline(zint,ig(:,10));
Sols11=spline(zint,ig(:,11));
for i=1:N
    Sols12{i}=spline(zint,ig(:,0*N+11+i));
    Sols13{i}=spline(zint,ig(:,1*N+11+i));
    Sols14{i}=spline(zint,ig(:,2*N+11+i));
    Sols15{i}=spline(zint,ig(:,3*N+11+i));
    Sols16{i}=spline(zint,ig(:,4*N+11+i));
    Sols17{i}=spline(zint,ig(:,5*N+11+i));
    if igs==0
        xiris{i}=spline(zint,xiri(:,i));
    else
        xiris{i}=spline(zint,xirig(:,i));
    end  
end
tstrt=tic;
for is0=1:ns0
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    s0=s00(is0)              % spray concentration at the wavecrest level
    pi3=3*D0*s0*zwt/alpv/us/r0^2
    pi5=3*K0*s0*zwt/alpT/us/r0^2/cpd/rhod0
    pi12=sigma*s0*g*zwt/us^2
    Solver;                  % Soving equations iteratively
    Sz1=ppval(Sols1,zint);   % ua
    Sz2=ppval(Sols2,zint);
    Sz3=ppval(Sols3,zint);   % ea
    Sz4=ppval(Sols4,zint);
    Sz5=ppval(Sols5,zint);   % epsl
    Sz6=ppval(Sols6,zint);
    Sz7=ppval(Sols7,zint);   % q
    Sz8=ppval(Sols8,zint);
    Sz9=ppval(Sols9,zint);   % Ta
    Sz10=ppval(Sols10,zint);
    Sz11=ppval(Sols11,zint); % P
    ka  =ppval(kas,zint);    % ka
    if is0==1 && ithq==0
        qinf=Sz7(end)
        thetainf=Sz9(end)/Sz11(end)^kd
    end
    Sz12=zeros(np,N);
    Sz13=Sz12; Sz14=Sz12; Sz15=Sz12; Sz16=Sz12; Sz17=Sz12; 
    xiri=Sz12; C0i=Sz12; C1i=Sz12; ai=Sz12; di=Sz12; ki=Sz12;
    for i=1:N
        Sz12(:,i)=ppval(Sols12{i},zint); % ui
        Sz13(:,i)=ppval(Sols13{i},zint);
        Sz14(:,i)=ppval(Sols14{i},zint); % lsi
        Sz15(:,i)=ppval(Sols15{i},zint);
        Sz16(:,i)=ppval(Sols16{i},zint); % Ti
        Sz17(:,i)=ppval(Sols17{i},zint);
        xiri(:,i)=ppval(xiris{i},zint);  % xir
        C0i(:,i) =ppval(C0is{i},zint);   % C0i
        C1i(:,i) =ppval(C1is{i},zint);   % C1i
        ai(:,i)  =ppval(ais{i},zint);    % ai
        di(:,i)  =ppval(dis{i},zint);    % di
        ki(:,i)  =ppval(kis{i},zint);    % ki
        end_
end
fprintf('Total time is %f seconds',toc(tstrt))
end
%============================ Main solver =================================
function Solver
global N np alp sig0t Cb1 Cb2
global pi0 pi9 pi17
global ri pii xiris C1is
global Sols1 Sols2 Sols3 Sols4 Sols5 Sols6 Sols7 Sols8 Sols9 Sols10
global Sols11 Sols12 Sols13 Sols14 Sols15 Sols16 Sols17
global kas xiris C0is C1is ais dis kis 
global options odeoptions tola tolit nitmax zint zplt iprst iter
eriter=1; iter=0;
fig1=figure('Position',[900  400 700 900]);
iprst=0; % set iprst=0 to compute persistent parameters
while eriter>tolit && iter<nitmax
    iter=iter+1
    if iter>1
        Sz=Szn;
    end
    %use old solution as an initial guess for a new iteration
    %   (see iguess routine below):
    solinit=bvpinit(zint,@iguess1);
    sol=bvp5c(@eqns,@bcs,solinit,odeoptions);
    Szn=deval(sol,zint);
    plt(N,zplt,zint,Szn,fig1)
    ua  =Szn(1,:)';
    ea  =Szn(3,:)';
    epsl=Szn(5,:)';
    q   =Szn(7,:)';
    Ta  =Szn(9,:)'; 
    P   =Szn(11,:)';
    ui  =Szn(12+0*N:11+1*N,:)';
    si  =exp(Szn(12+2*N:11+3*N,:))';
    dlsi=Szn(12+3*N:11+4*N,:)';
    % rhoa=P./Ta./(1+pi9*q);
    ka=alp^2*ea.^2./epsl;
%**** Solving nonlinear equations for coefficients & flagging complex
%**** solutions:
    xiri=.1*ones(np,N);
    for i=np:-1:1
        F=@(X)xir2(X,zint(i),ea(i),epsl(i),ua(i),ui(i,:),si(i,:),...
            dlsi(i,:));
        if i < np
            X=fsolve(F,xiri(i+1,:),options);
        else 
            X=fsolve(F,xiri(i,:),options);
        end
        for j=1:N
            if abs(imag(X(j)))>tola
                disp('Solver: Complex xir at'), zint(i)
                error('Stop Solver')
            end
        end
        xiri(i,:)=real(X);
    end
    for i=1:np
        for j=1:N
            if abs(xiri(i,j))<tola
                xiri(i,j)=max(xiri(:,j));
                %xiri(i,j)=1;
            end
        end
    end
    for j=1:N
        xiris{j}=spline(zint,xiri(:,j));
    end
%**** Comparing old and new solutions:
    if iter>1
       eriter=max(max(abs(Sz-Szn)))
    end
end
%spline coefficients for the solution:
Sols1=spline(zint,Szn(1,:)');
Sols2=spline(zint,Szn(2,:)');
Sols3=spline(zint,Szn(3,:)');
Sols4=spline(zint,Szn(4,:)');
Sols5=spline(zint,Szn(5,:)');
Sols6=spline(zint,Szn(6,:)');
Sols7=spline(zint,Szn(7,:)');
Sols8=spline(zint,Szn(8,:)');
Sols9=spline(zint,Szn(9,:)');
Sols10=spline(zint,Szn(10,:)');
Sols11=spline(zint,Szn(11,:)');
kas=spline(zint,ka);
wi=-diag(Hvs(pi0,zint))*(1./si)*diag(pii);
Cbi=Cb1+Cb2*(ui.^2+wi.^2-ui.*repmat(ua,1,N)).^2./(ui.^2+wi.^2)...
    ./((ui-repmat(ua,1,N)).^2+wi.^2);
Cbi(find((ui-repmat(ua,1,N)).^2+wi.^2<eps))=Cb1+Cb2;

C0i=1/sig0t./sqrt(1+Cbi.*xiri); uir=diag(sqrt(2*ea/3))*sqrt(xiri); 
ai=zeros(size(uir));
for i=1:np
    ai(i,:)=fai(ri,uir(i,:));
end
C1i=1./(1+2/3*pi17/alp^2*diag(epsl./ea)*ai./C0i);
di=diag(ka)*C0i; dai=C1i.*di;
ki=dai+pi17/3*diag(ea)*ai.*C1i;
for j=1:N
    Sols12{j}=spline(zint,Szn(11+0*N+j,:)');
    Sols13{j}=spline(zint,Szn(11+1*N+j,:)');
    Sols14{j}=spline(zint,Szn(11+2*N+j,:)');
    Sols15{j}=spline(zint,Szn(11+3*N+j,:)');
    Sols16{j}=spline(zint,Szn(11+4*N+j,:)');
    Sols17{j}=spline(zint,Szn(11+5*N+j,:)');
    C0is{j}  =spline(zint,C0i(:,j));
    C1is{j}  =spline(zint,C1i(:,j));
    ais{j}   =spline(zint,ai(:,j));
    dis{j}   =spline(zint,di(:,j));
    kis{j}   =spline(zint,ki(:,j));
end
if iter==nitmax
    disp('Solver: the maximum number of iterations is reached')
end
end
%============================ Equations ===================================
function dfdz=eqns(z,f)
%--------------------------------------------------------------------------
global iprst
persistent pr
if iprst==0
    pr=defpareqns;
    iprst=1;
end
N=pr{1}; alp=pr{2}; alpe=pr{3}; alpeps=pr{4}; sigma=pr{5}; s0=pr{6}; 
sig0t=pr{7}; Cb1=pr{8}; Cb2=pr{9}; pi0=pr{10}; pi1=pr{11}; pi2=pr{12}; 
pi3=pr{13}; pi5=pr{14}; pi6=pr{15}; pi7=pr{16}; pi8=pr{17}; pi9=pr{18};
pi12=pr{19}; pi13=pr{20}; pi17=pr{21}; l=pr{22}; U=pr{23}; nu=pr{24};
us=pr{25}; kd=pr{26}; q0=pr{27}; Ceps1=pr{28}; Ceps2=pr{29}; Ceps3=pr{30};
dr=pr{31}; ri=pr{32}; pii=pr{33}; xiris=pr{34}; r0=pr{35}; T0=pr{36};
D0=pr{37}; cpd=pr{38}; cpv=pr{39};

xiri=zeros(1,N); dfi=xiri;

ua=f(1); ea=f(3); epsl=f(5); q=f(7); Ta=f(9); P=f(11);
ui=f(12+0*N:11+1*N)'; lsi=f(12+2*N:11+3*N)'; si=exp(lsi); 
Ti=f(12+4*N:11+5*N)'; 
rhoa=P/Ta/(1+pi9*q); rho=rhoa+sigma*s0*sum(si); %rhod=rhoa*(1-q0*q);
rhoap=rhoa*(1+q*q0*(cpv/cpd-1));
ka=alp^2*ea^2/epsl;   %turbulent viscosity
uap=f(2)/ka/rhoa; eap=f(4)/alpe/ka/rhoa; epslp=f(6)/alpeps/ka/rhoa;
qp=f(8)/ka/rhoa; Tap=f(10)/ka/rhoap*P^kd-pi8*kd*rho*Ta/P;
for i=1:N
    xiri(i)=ppval(xiris{i},z);
end
wi=-pii./si*Hvs(pi0,z);

Cbi=Cb1+Cb2*(ui.^2+wi.^2-ui*ua).^2./(ui.^2+wi.^2)./((ui-ua).^2+wi.^2);
Cbi(find((ui-ua).^2+wi.^2<=eps))=Cb1+Cb2;


C0i=1/sig0t./sqrt(1+Cbi.*xiri); 
uir=sqrt(2*ea/3*xiri); ai=fai(ri,uir);
C1i=1./(1+2/3*pi17/alp^2*epsl/ea*ai./C0i);
di=C0i*ka; dai=C1i.*di;
lsip=(f(12+3*N:11+4*N)'-pi1*ai)./di;
wid=-dai.*lsip; wi=-pii./si*Hvs(pi0,z);
pa=ka*uap^2-pi12/rhoa*si*wid'...
    -pi13*ka*(pi9*qp/(1+pi9*q)+Tap/Ta+pi8*kd*rho/P);
qai=2*C1i*ea; qae=pi12/pi1*si./ai*(qai-2*ea+(pi1*wi-wid).*wid)';
paeps=Ceps1*epsl/ea*pa; daeps=Ceps2*epsl^2/ea; qaeps=Ceps3*epsl/ea*qae;

ki=dai+pi17/3*ai.*C1i*ea;
uip=f(12+1*N:11+2*N)'./ki;
Tip=(f(12+5*N:11+6*N)'-pi1*pii.*Ti./si*Hvs(pi0,z))./di; 
[pvsi,Di,Ki]=propTP(Ta,Ti,P,uir,N,nu,r0,T0,D0,ri,us);
dqi=rhoa*Di./ri.^2.*(q-pvsi./(P-pi9*pvsi));
if N>1
    dfi(1)=ri(2)/dr/6*(dqi(1)+si(2)/si(1)*dqi(2))-dqi(1);
    for i=2:N-1
        dfi(i)=-5/6*dqi(i)...
        +(ri(i+1)*si(i+1)*dqi(i+1)-ri(i)*si(i-1)*dqi(i-1))/si(i)/dr/6;            
    end
    dfi(N)=-ri(N)/dr/6*(si(N-1)/si(N)*dqi(N-1)+dqi(N))-dqi(N);
else
    dfi(1)=-dqi(1);
end
dfdz=[uap
      pi12/pi1*(si./ai)*(ua-ui)'
%
      eap
      rhoa*(epsl-pa)-qae
%
      epslp
      rhoa*(daeps-paeps)-qaeps
%
      qp
      pi3*dqi*si'
%
      Tap
      pi5*(Ki.*si./ri.^2)*(Ta-Ti)'
%
      -pi8*rho
%
      uip'
      ((ui-ua)/pi17./ai-ki.*uip.*lsip-pi1*pii.*uip./si*Hvs(pi0,z)...
                               +pi1*pii.*(ui-U)./si*delta(pi0,z))' 
%
      lsip'
      pi2*dfi'-f(12+3*N:11+4*N).*lsip'-pi1*(pii./si)'*delta(pi0,z)
%
      Tip'
      pi6*(Ki./ri.^2.*(Ti-Ta))'-(pi7*l*dqi+pi1*pii./si*delta(pi0,z))'-...
                               f(12+5*N:11+6*N).*lsip'
      ];
end

%========================= Boundary conditions ============================
function bc=bcs(fb,fi)
%Equations for
%f = [ua,ka*rhoa*ua',ea,ka*rhoa*ea',epsl,alpeps*ka*rhoa*epsl',...
%      q,ka*rhoa*q' ,Ta,ka*rhod*Theta',P,
%     ui,ki*ui',    ,ln(si),di*ln'(si)+pi1*ai,Ti,di*Ti'+pi1*H*pi*Ti/si]
%--------------------------------------------------------------------------
global N alp alpe alpeps sig0t Cb1 Cb2 pi0 pi1 pi8 pi9 pi17 kd
global tau ri xiris pii
global z1 zq zinf alpv alpT kp u1log q0 qs thetas T0 D0 nu r0 us
global is0 ithq thetainf qinf
% Values at the bottom edge of computational domain:
xirib=zeros(1,N); % dfib=xirib;
uab=fb(1); eab=fb(3); epslb=fb(5); qb=fb(7); Tab=fb(9); Pb=fb(11);
uib=fb(12+0*N:11+1*N)'; lsib=fb(12+2*N:11+3*N)'; sib=exp(lsib);
Tib=fb(12+4*N:11+5*N)';
kab=alp^2*eab^2/epslb; rhoab=Pb/Tab/(1+pi9*qb); 
for i=1:N
    xirib(i)=ppval(xiris{i},z1);
end
uirb=sqrt(2*eab/3*xirib); aib=fai(ri,uirb);
wib=-pii./sib*Hvs(pi0,z1);
Cbib=Cb1+Cb2*(uib.^2+wib.^2-uib*uab).^2./(uib.^2+wib.^2)...
    ./((uib-uab).^2+wib.^2);
Cbib(find((uib-uab).^2+wib.^2<=eps))=Cb1+Cb2;

C0ib=1/sig0t./sqrt(1+Cbib.*xirib); dib=C0ib*kab;
lsipb=(fb(12+3*N:11+4*N)'-pi1*aib)./dib; 
Tipb=(fb(12+5*N:11+6*N)'-pi1*pii.*Tib./sib*Hvs(pi0,z1))./dib;
[pvsib,Dib,~]=propTP(Tab,Tib,Pb,0,N,nu,r0,T0,D0,ri,us);
% [pvsib,Dib,Kib]=propTP(Tab,Tib,Pb,0,N,nu,r0,T0,D0,ri,us);
dqib=rhoab*Dib./ri.^2.*(qb-pvsib./(Pb-pi9*pvsib));

C1ib=1./(1+2/3*pi17/alp^2*epslb/eab*aib./C0ib);
daib=C1ib.*dib; kib=daib+pi17/3*aib.*C1ib*eab;
uipb=fb(12+1*N:11+2*N)'./kib;

% Values at infinity:
xirii=zeros(1,N);
uai=fi(1); eai=fi(3); epsli=fi(5); qi=fi(7); Tai=fi(9); Pi=fi(11);
uii=fi(12+0*N:11+1*N)'; lsii=fi(12+2*N:11+3*N)'; sii=exp(lsii);
kai=alp^2*eai^2/epsli; rhoai=Pi/Tai/(1+pi9*qi);
eapi=fi(4)/kai/rhoai/alpe; epslpi=fi(6)/kai/rhoai/alpeps;
for i=1:N
    xirii(i)=ppval(xiris{i},zinf);
end
wii=-pii./sii*Hvs(pi0,zinf);
Cbii=Cb1+Cb2*(uii.^2+wii.^2-uii*uai).^2./(uii.^2+wii.^2)...
    ./((uii-uai).^2+wii.^2);
Cbii(find(((uii-uai).^2+wii.^2)<=eps))=Cb1+Cb2;

C0ii=1/sig0t./sqrt(1+Cbii.*xirii); dii=C0ii*kai;
Tipi=fi(12+5*N:11+6*N)'./dii;

if is0==1 && ithq==0  % reference case
    bc=[
% ua:
        uab-u1log
        fi(2)-1
% ea:
        eab-1/alp
        eapi-pi8*(1-kd)/(tau^2-1)*tau^2/alp
% epsl:
        epslb-1/kp/z1
        epslpi+epsli/zinf
% q:
        qb-1/(1-pi9)+qs/alpv/kp/q0*log(z1/zq)
        fi(8)+qs/alpv/q0
% Ta:
        Tab/Pb^kd-1+thetas/alpT/kp/T0*log(z1/zq)
        fi(10)+thetas/alpT/T0
% P:
        Pb-1
% ui:
        %uib'-uab
        (uib-uab)/pi17./aib-pi1*pii.*uipb./sib
        uii'-uai
% si:
        lsipb'
        fi(12+3*N:11+4*N)
% Ti:
        Tipb'
        Tipi'
        ];
else
    qsb=-alpv*q0*fb(8); thetasb=-alpT*T0*fb(10);
    bc=[
% ua:
        uab-u1log
        fi(2)-1
% ea:
        eab-1/alp
        eapi-pi8*(1-kd)/(tau^2-1)*tau^2/alp
% epsl:
        epslb-1/kp/z1
        epslpi+epsli/zinf
% q:
        qb-1/(1-pi9)+qsb/alpv/kp/q0*log(z1/zq)
        qi-qinf
% Ta:
        Tab/Pb^kd-1+thetasb/alpT/kp/T0*log(z1/zq)
        Tai/Pi^kd-thetainf
% P:
        Pb-1
% ui:
        %uib'-uab
        (uib-uab)/pi17./aib-pi1*pii.*uipb./sib
        uii'-uai
% si:
        lsipb'
        fi(12+3*N:11+4*N)
% Ti:
        Tipb'
        Tipi'
        ];
end
end

%=========================== Initial guess ================================
function ic=iguess(z)
%**** Setting the initial guess: 
%--------------------------------------------------------------------------
global N kp alp alpeps z0 pi0 pi1
global s0i a0i
ic=zeros(11+6*N,1);
ic=[1/kp*log(z/z0)           % ua
    1                        % ka*ua'
    1/alp                    % ea
    0                        % alpe*ka*dea/dz
    1/kp/z                   % epsl
    -alpeps/z                % alpeps*ka*depsl/dz
    0.2                      % q
    0                        % ka*rhoa*dq/dz
    1.                       % Ta
    0                        % dTa/dz
    1                        % P
    1/kp*log(z/z0)*ones(N,1) % ui
    ones(N,1)                % ki*si*ui'
    -pi1/kp*log(s0i'*z)      % log(si)       
    pi1*(a0i'-1)             % k*dlog(si)/dz+pi1*ai
    ones(N,1)                % Ti
    pi1*a0i'*Hvs(pi0,z)];    % k*k*dTi/dz+pi1*ai*Hvs
end

%=========================== Initial guess 1 ==============================
function ic=iguess1(z)
%**** Interpolating the nitial guess:
%--------------------------------------------------------------------------
global N Sols1 Sols2 Sols3 Sols4 Sols5 Sols6 Sols7 Sols8 Sols9 Sols10
global Sols11 Sols12 Sols13 Sols14 Sols15 Sols16 Sols17
ic=zeros(11+6*N,1);
ic(1:11)=[ppval(Sols1,z); 
          ppval(Sols2,z); 
          ppval(Sols3,z); 
          ppval(Sols4,z); 
          ppval(Sols5,z); 
          ppval(Sols6,z);
          ppval(Sols7,z);
          ppval(Sols8,z);
          ppval(Sols9,z);
          ppval(Sols10,z);
          ppval(Sols11,z);
          ];
for i=1:N
    ic(0*N+11+i)=ppval(Sols12{i},z);
    ic(1*N+11+i)=ppval(Sols13{i},z);
    ic(2*N+11+i)=ppval(Sols14{i},z);
    ic(3*N+11+i)=ppval(Sols15{i},z);
    ic(4*N+11+i)=ppval(Sols16{i},z);
    ic(5*N+11+i)=ppval(Sols17{i},z);
end
end

%======================= Local terminal velocity ==========================
function ai=fai(ri,uir)
global sigma g r0 a0 us nu
Ru=2*r0*us/nu*ri.*uir;
cdi=(3.17e8+6.69e7*Ru+4.47e5*Ru.^2)./(1.37e7*Ru+9.86e5.*Ru.^2-Ru.^3);
ai=8*sigma*g*r0/3/us/a0*ri./cdi./uir;
end

%========================= Physical properties ============================
function [pvsi,Di,Ki]=propTP(Ta,Ti,P,uir,N,nu,r0,T0,D0,ri,us)
fvi=ones(size(uir));
%******* Thermal conductivity
K=(1-4.69e-3*T0*Ta+1.08e-5*T0^2*Ta^2)/(1-4.69e-3*T0+1.08e-5*T0^2);
%******* Vapour diffusion coefficient
D=(1-7.00e-3*T0*Ta+1.07e-5*T0^2*Ta^2)/(1-7.00e-3*T0+1.07e-5*T0^2)/P;
if uir>0
%******* Ventilation coefficients
    Rei=2*us*uir.*ri*r0/nu; Sc=nu/D0/D; chi=Sc^(1/3)*sqrt(Rei);
    fvi=zeros(1,N);
    for j=1:N
        if chi(j)<=1.4
            fvi(j)=1+.108*chi(j)^2;
        else
            fvi(j)=0.78048+.308*chi(j);
        end
    end 
end
Di=fvi*D; Ki=fvi*K;
%******* Saturated vapour pressure
pvsi=Ti.^(-5.09).*exp(6808*(Ti-1)/T0./Ti);
end

%==========================================================================
function F=termv(a,sigma,nu,rav,g)
%terminal velocity of a droplet of radius rav in still air
Re=2*a*rav/nu;
Cd=3808*(1617933/2030+178861/1063*Re+1219/1084*Re^2)/...
    (681*Re*(77531/422+13529/976*Re-Re^2/71154));
F=3/8*Cd/sigma*a^2/rav-g;
end

%==========================================================================
function dlt=delta(pi0,z)
dlt=0;
if abs(z-1)<=1/pi0
   dlt=pi0/2*(1+cos(pi*pi0*(z-1)));
end
end

%==========================================================================
function H=Hvs(pi0,z)
H=1.;
if abs(z-1)<=1/pi0
   H=(1-pi0*(z-1)-sin(pi*pi0*(z-1))/pi)/2;
end
if z-1>1/pi0
   H=0;
end
end

%==========================================================================
function dFmdr=Andreas98(ri,us,kp,u1log)
%--------------------------------------------------------------------------
p=[0.065 0.49 0 -1000*us^2];            % solving (2.4) and (2.5) for U10
buf=roots(p);
for i=1:3
    if imag(buf(i))==0
        U10=real(buf(i));
    end
end
U10=us*(u1log+log(2)/kp);                % logarithmic U10
if us<sqrt(0.12)*11/10
    CDN10=1.205e-3;                      % (2.5a)
elseif us>=sqrt(0.12)*11/10
    CDN10=(0.49+0.065*U10)*1.e-3;        % (2.5b)
else
    error('Andreas98 error:', 'U10 is out of range');
end
U14=U10*(1+sqrt(CDN10)*log(1.4)/kp);     % (A5)
A1=10^(0.0676*U14+2.43);                 % (A2a)
A2=10^(0.959*sqrt(U14)-1.476);           % (A2b)
r1=2.1; r2=9.2; f1=3.1; f2=3.3;          
t1=f1*(log(10/r1))^2;
t2=f2*(log(10/r2))^2;
CA1=10*(A1*exp(-t1)+A2*exp(-t2));        % (A1), explanation after (3.5)
CA2=CA1*37.5^1.8;
CA3=CA2*100^5.2;
r80=0.518*(10^6*ri).^(0.976);            % (3.4), r_0=10^6*ri is in microns
for i=1:length(ri)
    if r80(i)<10
        t1=f1*(log(r80(i)/r1))^2;
        t2=f2*(log(r80(i)/r2))^2;
        dFsdr(i)=A1*exp(-t1)+A2*exp(-t2);
    elseif r80(i)>=10 && r80(i)<37.5     % (3.5)
        dFsdr(i)=CA1/r80(i);
    elseif r80(i)>=37.5 && r80(i)<100
        dFsdr(i)=CA2*r80(i)^(-2.8);
    elseif r80(i)>=100 && r80(i)<=2000
        dFsdr(i)=CA3*r80(i)^(-8);
    else
        error('Andreas98 error: r80 is out of range');
    end
end
drsdr=0.506*(10^6*ri).^(-0.024);         % (3.7), r_0=10^6*ri is in microns
dFmdr=3.5*(4*pi/3)*ri.^3.*dFsdr.*drsdr;  % (3.8)
end

%==========================================================================
function dFmdr=OrtizSuslow16(ri)
%--------------------------------------------------------------------------
p=[-2.1009 7.3857 -3.8341];
R=6+log10(ri);
dFmdr=4/3*pi*ri.^3.*10.^(p(1)*R.^2+p(2)*R+p(3));
end

%==========================================================================
function z0t=roughl(us,iz0,g,kp)
zr=10; %m, reference altitude
if iz0==0
    alp_C=0.015;
else
    f=@(x) (0.49+0.065*us*x)*x^2-1000;
    options=optimset('TolFun',1.e-10);
    x=fsolve(f,10,options);
    alp_C=g*zr*exp(-x*kp)/us^2; %Charnock constant
end
z0t=alp_C*us^2/g; %roughness length
end

%==========================================================================
function F=xir2(xiri,z,ea,epsl,ua,ui,si,dlsi)
%**** Equation for $xi_r^2$
global sig0t alp pi0 pi1 pi17 Cb1 Cb2 ri pii
wi=-pii./si*Hvs(pi0,z);
Cbi=Cb1+Cb2*(ui.^2+wi.^2-ui*ua).^2./(ui.^2+wi.^2)./((ui-ua).^2+wi.^2);
Cbi(find((ui-ua).^2+wi.^2<=eps))=Cb1+Cb2;

uir=sqrt(2*ea/3*xiri); C0i=1/sig0t./sqrt(1+Cbi.*xiri); ai=fai(ri,uir);
C1i=1./(1+2*pi17/3/alp^2*epsl/ea*ai./C0i);
dlsidz=C1i.*(dlsi-pi1*ai);
F=2*ea*xiri/3-(ui-ua).^2-(pi1*wi+dlsidz).^2-2*ea*(1-C1i);
end

%==========================================================================
function plt(N,zplt,zint,Szn,fig1)
global pi0 pi9 pi17 alp kp kd z0 sig0t Cb1 Cb2
global options xiris ri pii tola
np=length(zint);
xiri=zeros(np,N); ai=xiri;
fsz=16; lw=2;
ua  =Szn(1,:)';
ea  =Szn(3,:)';
epsl=Szn(5,:)';
q   =Szn(7,:)';
Ta  =Szn(9,:)'; 
P   =Szn(11,:)';
ui  =Szn(12+0*N:11+1*N,:)';
si  =exp(Szn(12+2*N:11+3*N,:)');
dlsi=Szn(12+3*N:11+4*N,:)';
Ti  =Szn(12+4*N:11+5*N,:)';
rhoa=P./Ta./(1+pi9*q);
ka=alp^2*ea.^2./epsl;
tet=Ta./P.^kd;
tetv=(1+pi9*q).*tet;

for j=1:N
    xiri(:,j)=ppval(xiris{j},zint);
end
uir=diag(sqrt(2*ea/3))*sqrt(xiri);
for i=1:np
    ai(i,:)=fai(ri,uir(i,:));
end
wi=-diag(Hvs(pi0,zint))*(1./si)*diag(pii);
Cbi=Cb1+Cb2*(ui.^2+wi.^2-ui.*repmat(ua,1,N)).^2./(ui.^2+wi.^2)...
    ./((ui-repmat(ua,1,N)).^2+wi.^2);
Cbi(find((ui-repmat(ua,1,N)).^2+wi.^2<-eps))=Cb1+Cb2;
C0i=1/sig0t./sqrt(1+Cbi.*xiri);
C1i=1./(1+2/3*pi17/alp^2*diag(epsl./ea)*ai./C0i);
di=diag(ka)*C0i; dai=C1i.*di; ki=dai+pi17/3*diag(ea)*ai.*C1i;
xirin=.1*ones(np,N);
for i=np:-1:1
    F=@(X)xir2(X,zint(i),ea(i),epsl(i),ua(i),ui(i,:),si(i,:),dlsi(i,:));
    if i < np
        xirin(i,:)=fsolve(F,xirin(i+1,:),options);
    else
        xirin(i,:)=fsolve(F,xirin(i,:),options);
    end
    for j=1:N
        if abs(imag(xirin(i,j)))>tola
            disp('Solver: Complex coefficients at'), zint(i)
            error('Stop plt')
        end
    end
end
uirn=diag(sqrt(2*ea/3))*sqrt(xirin);
for i=1:np
    ain(i,:)=fai(ri,uirn(i,:));
end
C0in=1/sig0t./sqrt(1+Cbi.*xirin); 
C1in=1./(1+2/3*pi17/alp^2*diag(epsl./ea)*ain./C0in);
din=diag(ka)*C0in; dain=C1in.*din; kin=dain+pi17/3*diag(ea)*ain.*C1in;

%--------------------------------------------------------------------------
ip=find(zint<=zplt+.1);
figure(fig1)
subplot(4,2,1)
plot(zint(ip),ua(ip),'LineWidth',lw), hold on, box on
ulog=1/kp*log(zint/z0);
plot(zint(ip),ulog(ip),'--g','Linewidth',lw)
xlabel('$z$','Interpreter','latex','FontSize',fsz);
ylabel('$u_a$','Interpreter','latex','FontSize',fsz);
axis('tight')
hold off

subplot(4,2,2)
plot(zint(ip),ea(ip),'LineWidth',lw), hold on, box on
xlabel('$z$','Interpreter','latex','FontSize',fsz);
ylabel('$e_a$','Interpreter','latex','FontSize',fsz);
axis('tight')
hold off

subplot(4,2,3)
plot(zint(ip),epsl(ip),'LineWidth',lw), hold on, box on
xlabel('$z$','Interpreter','latex','FontSize',fsz);
ylabel('$\epsilon_a$','Interpreter','latex','FontSize',fsz);
axis('tight')
hold off

subplot(4,2,4)
plot(zint(ip),q(ip),'LineWidth',lw), hold on, box on
xlabel('$z$','Interpreter','latex','FontSize',fsz);
ylabel('$q$','Interpreter','latex','FontSize',fsz);
axis('tight')
hold off

subplot(4,2,5)
plot(zint(ip),Ta(ip),'-','LineWidth',lw), hold on, box on
plot(zint(ip),tetv(ip),'--','LineWidth',lw)
plot(zint(ip),tet(ip),'-.','LineWidth',lw)
xlabel('$z$','Interpreter','latex','FontSize',fsz);
ylabel('$T_a,\ \theta,\ \theta_v$','Interpreter','latex','FontSize',fsz);
axis('tight')
hold off

subplot(4,2,6)
plot(zint(ip),P(ip),'LineWidth',lw), hold on, box on
xlabel('$z$','Interpreter','latex','FontSize',fsz);
ylabel('$P$','Interpreter','latex','FontSize',fsz);
axis('tight')
hold off

subplot(4,2,7)
plot(zint(ip),rhoa(ip),'LineWidth',lw), hold on, box on
xlabel('$z$','Interpreter','latex','FontSize',fsz);
ylabel('$\rho_a$','Interpreter','latex','FontSize',fsz);
axis('tight')
hold off

subplot(4,2,8)
plot(zint(ip),ka(ip),'LineWidth',lw), hold on, box on
plot(zint(ip),kp*zint(ip),'g--','Linewidth',lw)
xlabel('$z$','Interpreter','latex','FontSize',fsz);
ylabel('$k_a$','Interpreter','latex','FontSize',fsz);
axis('tight')
hold off

pause(1)
end

%===================Define values of parameters for eqns ==================
function pr=defpareqns
global N alp alpe alpeps sigma s0 r0 sig0t Cb1 Cb2
global pi0 pi1 pi2 pi3 pi5 pi6 pi7 pi8 pi9 pi12 pi13 pi17 l U cpd cpv
global nu us kd q0 Ceps1 Ceps2 Ceps3
global dr ri pii xiris T0 D0
pr={N,alp,alpe,alpeps,sigma,s0,sig0t,Cb1,Cb2,pi0,pi1,pi2,pi3,pi5,pi6,...
    pi7,pi8,pi9,pi12,pi13,pi17,l,U,nu,us,kd,q0,Ceps1,Ceps2,Ceps3,dr,ri,...
    pii,xiris,r0,T0,D0,cpd,cpv};
end
%==========================================================================