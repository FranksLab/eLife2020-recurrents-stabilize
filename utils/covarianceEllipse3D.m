function covarianceEllipse3D(c,P,v)
% covarianceEllipse3D plots a Gaussian as an uncertainty ellipse
% Based on Maybeck Vol 1, page 366
% k=2.296 corresponds to 1 std, 68.26% of all probability
% k=11.82 corresponds to 3 std, 99.74% of all probability
%

%% this code i got from github
% Modified from http://www.mathworks.com/matlabcentral/newsreader/view_thread/42966
% 
% [e,s] = svd(P);
% k = 2.296; 
% radii = k*sqrt(diag(s));
% 

%% this code is from the thread linked by the above code
% For N standard deviations spread of data, the radii of the eliipsoid will
% be given by N*SQRT(eigenvalues).
M = c;
[U,L] = eig(P);

N = 1; % choose your own N
radii = N*sqrt(diag(L));

% generate data for "unrotated" ellipsoid
[xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));

% rotate data with orientation matrix U and center M
a = kron(U(:,1),xc); b = kron(U(:,2),yc); c = kron(U(:,3),zc);
data = a+b+c;  n = size(data,2);
x = data(1:n,:)+M(1); 
y = data(n+1:2*n,:)+M(2); 
z = data(2*n+1:end,:)+M(3);


% now plot the rotated ellipse
% sc = surf(x,y,z,abs(xc));
%         surf(x,y,z,v*ones(size(x)),'FaceLighting','gouraud'); %shading interp; 
        p = surf(x,y,z,v*ones(size(x))); %shading interp; 


shading interp
alpha(0.36)
axis equal