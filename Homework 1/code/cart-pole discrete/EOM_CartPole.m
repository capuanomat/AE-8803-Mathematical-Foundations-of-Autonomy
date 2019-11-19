 function [dXdt, xdot_x, xdot_u] = EOM_CartPole(x, u)
 
EPS = 1e-5;

global E

mp    = E.mp;
mc    = E.mc;
muc   = E.muc;
mup   = E.mup;
l     = E.l;
g     = E.g;

V     = x(2,:);
TH    = x(3,:);
OM    = x(4,:);

Q     = 1./(8*mc + (5 - 3*cos(2*TH))*mp);
A     = sin(u/3)/1.5 + u - muc*V - mp*sin(TH)*l.*OM.^2;
B     = g*sin(TH)*mp-mup*OM/l;

AccX  = (8*A + 6*cos(TH).*B).*Q;
Z1    = 6*(mp*cos(TH).*A + (mc+mp)*B);
AccTH = Z1.*(Q/(l*mp));

dXdt  = [V; AccX; OM; AccTH];

%----------- compute xdot_x using finite differences ------------
if nargout>1,
    
    nx = size(x,1);
    nu = size(u,1);
    
    x1 = repmat(x, [1,nx]) + eye(nx)*EPS;
    x2 = repmat(x, [1,nx]) - eye(nx)*EPS;
    uu = repmat(u,[1,nx]);
    
    f1 = EOM_CartPole(x1,uu);
    f2 = EOM_CartPole(x2,uu);
    
    xdot_x = (f1-f2)/2/EPS;
    
    u1 = repmat(u, [1,nu]) + eye(nu)*EPS;
    u2 = repmat(u, [1,nu]) - eye(nu)*EPS;
    xx = repmat(x,[1,nu]);
    
    g1 = EOM_CartPole(xx,u1);
    g2 = EOM_CartPole(xx,u2);
    
    xdot_u = (g1-g2)/2/EPS;
  
end

end

