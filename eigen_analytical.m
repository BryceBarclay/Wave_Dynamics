function [PHI,lam] = eigen_analytical(numeig,x,y,z,L,W,H)
% eigenvalues and eigenfunctions: for domain [0,L]X[0,W]X[0,H]

nx = length(x);
ny = length(y);
nz = length(z);
[ll,mm,nn] = ndgrid((pi/L)*(1:numeig),(pi/W)*(1:numeig),(pi/H)*(1:numeig));  % 1D eigenvalues
[XX,YY,ZZ] = ndgrid(x,y,z);
ll = reshape(ll,[numeig^3,1]); mm = reshape(mm,[numeig^3,1]); nn = reshape(nn,[numeig^3,1]);
lam = ll.^2 + mm.^2 + nn.^2;                    % 3D eigenvalues
lam = reshape(lam,[numeig^3,1]);

[lam,I] = sort(lam);
ll = ll(I);
mm = mm(I);
nn = nn(I);

% Compute eigenfunctions as orthonormal basis
PHI = zeros(nx,ny,nz,numeig^3);
for k = 1:numeig^3
    PHI(:,:,:,k) = (8/(L*W*H))^0.5 * sin(ll(k)*XX).*sin(mm(k)*YY).*sin(nn(k)*ZZ);
end

end