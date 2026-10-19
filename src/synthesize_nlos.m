function [Y,info]=synthesize_nlos(cfg,s,seed)
% Narrowband adaptation of BMT's unit-hop template and received-total kappa.
% Cascade contains LoS*NLoS + NLoS*LoS + NLoS*NLoS. The separate single
% bounce follows user -> fixed scatterer -> BS and bypasses every RIS.
Y=zeros(size(s.Y_los)); info=struct('input_kappa_db',Inf,'projected_kappa_db',Inf, ...
    'raw_cascade_power',0,'raw_direct_power',0,'mixed_cascade_power',0, ...
    'mixed_direct_power',0,'cross_power',0,'mixed_power',0,'scale',0);
if strcmp(cfg.nlos_model,'los') || isinf(cfg.kappa_db), return; end
assert(ismember(cfg.nlos_model,{'cascade10','cascade20','direct20','mixed20'}));
L=20; if strcmp(cfg.nlos_model,'cascade10'), L=10; end
% Always draw a 20-path master template, including the disabled branch.
% Prefix selection and branch switches cannot shift any other random draw.
r=RandStream('mt19937ar','Seed',seed+8000000);
t.bs_dir=zeros(3,20,cfg.RIS_m); t.br_dir=t.bs_dir;
t.ru_dir=zeros(3,20,cfg.User_k,cfg.RIS_m);
t.br_gain=(randn(r,20,cfg.RIS_m)+1j*randn(r,20,cfg.RIS_m))/sqrt(2);
t.ru_gain=(randn(r,20,cfg.User_k,cfg.RIS_m)+1j*randn(r,20,cfg.User_k,cfg.RIS_m))/sqrt(2);
for m=1:cfg.RIS_m
    t.bs_dir(:,:,m)=directions(r); t.br_dir(:,:,m)=directions(r);
    for k=1:cfg.User_k, t.ru_dir(:,:,k,m)=directions(r); end
end
t.delay=1e-7*rand(r,20,cfg.User_k);
t.direct_gain=(randn(r,20,cfg.User_k)+1j*randn(r,20,cfg.User_k))/sqrt(2);
t.direct_dir=zeros(3,20,cfg.User_k);
for k=1:cfg.User_k, t.direct_dir(:,:,k)=directions(r); end
C=zeros(size(Y)); D=zeros(size(Y));
if ~strcmp(cfg.nlos_model,'direct20')
    for m=1:cfg.RIS_m
        ar=array2D(s.phi_RIS_x(m),s.phi_RIS_y(m),cfg.RIS_size_x,cfg.RIS_size_y);
        brLos=s.lambda/(4*pi*s.d_B(m))*exp(-1j*2*pi*s.d_B(m)/s.lambda)*(s.A_BS(:,m)*ar');
        ab=array((pi*s.BS_dir.'*t.bs_dir(:,1:L,m)).',cfg.BS_size);
        arN=array2D((pi*s.RIS_dir_x(:,m).'*t.br_dir(:,1:L,m)).', ...
            (pi*s.RIS_dir_y(:,m).'*t.br_dir(:,1:L,m)).',cfg.RIS_size_x,cfg.RIS_size_y);
        brN=(ab.*t.br_gain(1:L,m).')*arN';
        brN=brN*(norm(brLos,'fro')/norm(brN,'fro'));
        physicalOmega=s.Omega(:,:,m).*ar;
        for k=1:cfg.User_k
            ruLos=s.lambda/(4*pi*s.d_U(k,m))*exp(-1j*2*pi*s.d_U(k,m)/s.lambda) ...
                *array2D(s.th_RIS_x(k,m),s.th_RIS_y(k,m),cfg.RIS_size_x,cfg.RIS_size_y);
            aN=array2D((pi*s.RIS_dir_x(:,m).'*t.ru_dir(:,1:L,k,m)).', ...
                (pi*s.RIS_dir_y(:,m).'*t.ru_dir(:,1:L,k,m)).',cfg.RIS_size_x,cfg.RIS_size_y);
            ruN=aN*t.ru_gain(1:L,k,m); ruN=ruN*(norm(ruLos)/norm(ruN));
            % Matrix products include all L + L + L^2 diffuse combinations.
            C=C+sqrt(s.Tx_power)*(brLos*(physicalOmega.*ruN) ...
                +brN*(physicalOmega.*(ruLos+ruN))).*s.X(:,k).';
        end
    end
end
points=zeros(3,20,cfg.User_k);
for k=1:cfg.User_k
    offset=s.User_pos(:,k)-s.BS_pos; distance=norm(offset);
    excess=3e8*t.delay(:,k).'; total=distance+excess;
    radius=excess.*(2*distance+excess)./(2*(total-offset.'*t.direct_dir(:,:,k)));
    points(:,:,k)=s.BS_pos+t.direct_dir(:,:,k).*radius;
    total=vecnorm(points(:,:,k)-s.BS_pos)+vecnorm(points(:,:,k)-s.User_pos(:,k));
    gain=t.direct_gain(:,k)./total.'.*exp(-1j*2*pi*total.'/s.lambda);
    a=array((pi*s.BS_dir.'*t.direct_dir(:,:,k)).',cfg.BS_size);
    D=D+sqrt(s.Tx_power)*(a*gain)*s.X(:,k).';
end
pc=mean_power(C); pd=mean_power(D); Cm=zeros(size(C)); Dm=zeros(size(D));
switch cfg.nlos_model
    case 'mixed20', Cm=C/sqrt(2); Dm=D*sqrt(pc/(2*pd));
    case 'direct20', Dm=D;
    otherwise, Cm=C;
end
mixed=Cm+Dm; pm=mean_power(mixed); assert(isfinite(pm)&&pm>realmin,'Degenerate NLoS realization.');
scale=sqrt(mean_power(s.Y_los)/(10^(cfg.kappa_db/10)*pm)); Y=scale*mixed;
projection=(s.A_BS'*s.A_BS)\s.A_BS';
info.input_kappa_db=10*log10(mean_power(s.Y_los)/mean_power(Y));
info.projected_kappa_db=10*log10(mean_power(projection*s.Y_los)/mean_power(projection*Y));
info.raw_cascade_power=pc; info.raw_direct_power=pd;
info.mixed_cascade_power=mean_power(Cm); info.mixed_direct_power=mean_power(Dm);
info.cross_power=2*real(mean(conj(Cm(:)).*Dm(:))); info.mixed_power=pm; info.scale=scale;
% Diagnostics stay in this scene only; formal records store their scalars.
info.raw_direct=D; info.scatterer_position=points; info.excess_delay=t.delay; info.template=t;
end
function d=directions(r)
z=2*rand(r,1,20)-1; az=2*pi*rand(r,1,20); radius=sqrt(max(0,1-z.^2));
d=[radius.*cos(az);radius.*sin(az);z];
end
function p=mean_power(x), p=mean(abs(x(:)).^2); end
