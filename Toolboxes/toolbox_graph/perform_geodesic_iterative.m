function [U,err,Usvg] = perform_geodesic_iterative(vertex, faces, W, I, options)

options.null = 0;
verb = getoptions(options, 'verb', 1);
svg_rate = getoptions(options, 'svg_rate', 0);  % Set to 0 to disable Usvg

% ✅ Fix for infinite loop warning
if isfield(options, 'niter') && isnumeric(options.niter) && isfinite(options.niter)
    niter = options.niter;
else
    niter = 200;
end

dotp = @(u,v)sum(u.*v,1);
R = @(u)reshape(u, [1 1 length(u)]);
Inv1 = @(M,d)[M(2,2,:)./d -M(1,2,:)./d; -M(2,1,:)./d M(1,1,:)./d];
Inv  = @(M)Inv1(M, M(1,1,:).*M(2,2,:) - M(1,2,:).*M(2,1,:));
Mult = @(M,u)[M(1,1,:).*u(1,1,:) + M(1,2,:).*u(2,1,:);  M(2,1,:).*u(1,1,:) + M(2,2,:).*u(2,1,:)];

n = size(vertex,2);

i = [faces(1,:) faces(2,:) faces(3,:)];
j = [faces(2,:) faces(3,:) faces(1,:)];
k = [faces(3,:) faces(1,:) faces(2,:)];

err = [];
U = getoptions(options, 'U', []);
if isempty(U)
    U = zeros(n,1);
end

% Disable Usvg unless svg_rate > 0
if svg_rate > 0
    Usvg = zeros(n, 0);
else
    Usvg = [];
end

x  = vertex(:,i);
x1 = vertex(:,j) - x;
x2 = vertex(:,k) - x;

C = [R(dotp(x1,x1)) R(dotp(x1,x2)); ...
    R(dotp(x2,x1)) R(dotp(x2,x2))];
S = Inv(C);
a = sum(sum(S));
w = R(W(i));

L1 = sqrt(dotp(x1,x1)); L1 = L1(:).*w(:); 
L2 = sqrt(dotp(x2,x2)); L2 = L2(:).*w(:); 

for it=1:niter
    if verb
        progressbar(it,niter);
    end
    uj = U(j);
    uk = U(k);
    u = [R(uj); R(uk)];
    
    b = dotp(sum(S,2), u);
    c = dotp(Mult(S,u), u) - w.^2;
    delta = max(b.^2 - a.*c, 0);
    d = (b + sqrt(delta))./a;
    
    alpha = Mult(S, u - repmat(d, 2, 1));
    J = find(alpha(1,1,:)>0 | alpha(2,1,:)>0);
    
    d1 = L1 + uj(:);
    d2 = L2 + uk(:);
    d = d(:);
    d(J) = min(d1(J), d2(J));
    
    U1 = accumarray(i', d, [n 1], @min);  U1(U1==0) = Inf;
    U1(I) = 0;
    
    err(end+1) = norm(U-U1, 'fro');
    if err(end)==0
        break;
    end
    
    U = U1;
    
    % OPTIONAL: Save snapshots only if svg_rate > 0
    if svg_rate > 0 && mod(it, svg_rate)==0 && nargout > 2
        Usvg(:,end+1) = U;
    end
end
