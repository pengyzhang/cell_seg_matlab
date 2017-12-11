%% @Pengyue Zhang
%% Function for level set optimization
function [u, allU, cu, cb, transform] = lse(u0, allU, Img, g, R, peak, transform, lambdaU, lambdaB, mu, xi, omega, nu, timestep, epsilon, iter_lse)

u=u0;
[hatx,alignedDict] = updateSR(u, Img, R, transform);

[cu, cb] = updatef(u, allU, Img, epsilon);

[u, allU, transform] = updateLSF(u, allU, Img, cu, cb, g, R, alignedDict, hatx, peak, transform, lambdaU, lambdaB, mu, xi, omega, nu, timestep, epsilon, iter_lse);


% update sparse representation
function [hatx,alignedDict] = updateSR(u, Img, dict, transform)
nPhi = size(u,2);
[nRow, nCol] = size(Img);
hatx = [];
for iPhi = 1:nPhi
    
    t = transform{iPhi};
    [dim,nSample] = size(dict);
    [nRow,nCol] = size(Img);
    [dmDict randomMat]= distance_map(Img, dict, t);
    vecDict = [];
    for iSample = 1:nSample
        dmVec = dmDict{iSample};
        dmVec = reshape(dmVec,nRow*nCol,1);
        dmVec = randomMat*dmVec;
        dmVec = dmVec/norm(dmVec);
        vecDict = [vecDict dmVec];
    end
    
    y = reshape(u{iPhi},nRow*nCol,1);
    y = randomMat*y;
    y = y/norm(y);
    
    
    [tempx,e] = sparse_solver(y,vecDict);
    hatx = [hatx tempx];
    alignedDict{iPhi} = dmDict;
end

% update cu, cb
function [cu, cb] = updatef(u, allU, Img, epsilon)

nPhi = size(u,2);

for i = 1:nPhi
    Hu=Heaviside(u{i},epsilon);
    NuMat = Img.*Hu;
    Nu(i) = sum(NuMat(:));
    Du(i) = sum(Hu(:));
    cu(i) = Nu(i)/(Du(i)+eps);
end
Hb = 1-Heaviside(allU,epsilon);
NbMat = Img.*Hb;
Nb = sum(NbMat(:));
Db = sum(Hb(:));
cb = Nb/(Db+eps);

% update level set function
function [u, allU, transform] = updateLSF(u, allU, Img, cu, cb, g, R, alignedDict, hatx, peak, transform, lambdaU, lambdaB, mu, xi, omega, nu, timestep, epsilon, iter_lse)%u = updateLSF(Img, u0, C, KONE_Img, KB1, KB2, mu, nu, timestep, epsilon, iter_lse)

nPhi = size(u,2);

for iIter=1:iter_lse
    for iPhi = 1:nPhi
        u{iPhi}=NeumannBoundCond(u{iPhi});
        B = ones(size(Img))-(Heaviside(allU,epsilon)-Heaviside(u{iPhi},epsilon));
        DiracU=Dirac(u{iPhi},epsilon);
        e1 = (Img-ones(size(Img))*cu(iPhi)).^2;
        e2 = (Img-ones(size(Img))*cb).^2.*B;
        e1Term = lambdaU*e1;
        e2Term = lambdaB*e2;
        ImageTerm=-DiracU.*(e1Term-e2Term); % e1.*g e2.*P
        
        K=div_norm(u{iPhi});    % div( (\nabla u) / |u| )
        [gx,gy] = gradient(g);
        [ux,uy] = gradient(u{iPhi});
        normDu = sqrt(ux.^2+uy.^2+1e-10);
        Nx = ux./(normDu+eps);
        Ny = uy./(normDu+eps);
        lengthTerm = mu*DiracU.*(gx.*Nx+gy.*Ny+g.*K);
        %penalizeTerm=beta*(4*del2(u)-K);
         
        penalizeTerm = xi*distReg_p2(u{iPhi});
        exclusiveTerm = -omega*DiracU.*(ones(size(Img))-B);
        
        y = u{iPhi};
        y = y/norm(y,'fro');
        
        dmDict = alignedDict{iPhi};
        nSample = size(dmDict,2);
        sparseCoef = hatx(:,iPhi);
        sparseTerm = -2*nu*(y-sparse_recon(dmDict, sparseCoef));
            
        u{iPhi}=u{iPhi}+timestep*(ImageTerm+lengthTerm+penalizeTerm+exclusiveTerm+sparseTerm);
        peakX = peak(iPhi,1);
        peakY = peak(iPhi,2);
        u{iPhi} = post_process(u{iPhi},peakX, peakY);
        
        theta = acos(transform{iPhi}.T(1,1));
        s = transform{iPhi}.b;
        xIdx = transform{iPhi}.c(1,1);
        yIdx = transform{iPhi}.c(1,2);
        
        
        thetaGradient = zeros(size(Img));
        thetaGradient1 = zeros(size(Img));
        thetaGradient2 = zeros(size(Img));
        for iX1 = 1:size(Img,2)
            for iX2 = 1:size(Img,1)
                thetaGradient1(iX2,iX1) = s*(-iX1*sin(theta)+iX2*cos(theta));
                thetaGradient2(iX2,iX1) = s*(-iX1*cos(theta)-iX2*sin(theta));
            end
        end
        for iSample = 1:nSample
            [phiX, phiY] = gradient(dmDict{iSample});
            
            thetaGradient = thetaGradient+(phiX.*thetaGradient1+phiY.*thetaGradient2)*sparseCoef(iSample);
        end
        theta2D = 2*nu*(y-sparse_recon(dmDict, sparseCoef)).*(thetaGradient);
        thetaSum = sum(theta2D(:));
        theta = theta+timestep*thetaSum;
        transform{iPhi}.T = [cos(theta),sin(theta);-sin(theta),cos(theta)];
        
        xGradient = zeros(size(Img));
        yGradient = zeros(size(Img));
        for iSample = 1:nSample
            [phiX, phiY] = gradient(dmDict{iSample});
            xGradient = xGradient+phiX*sparseCoef(iSample);
            yGradient = yGradient+phiY*sparseCoef(iSample);
        end
        x2D = 2*nu*(y-sparse_recon(dmDict, sparseCoef)).*xGradient;
        xSum = sum(x2D(:));
        xIdx = xIdx+timestep*xSum;
        transform{iPhi}.c(:,1) = ones(size(R,1)/2,1)*xIdx;
    
        y2D = 2*nu*(y-sparse_recon(dmDict, sparseCoef)).*yGradient;
        ySum = sum(y2D(:));
        yIdx = yIdx+timestep*ySum;
        transform{iPhi}.c(:,2) = ones(size(R,1)/2,1)*yIdx;

    end
    c0 = 5;
    allU = -c0*ones(size(Img));
    for iPhi = 1:nPhi
        utemp = u{iPhi};
        allU(utemp>0) = c0;

    end
end

function f = distReg_p2(phi)
% compute the distance regularization term with the double-well potential p2 in eqaution (16)
[phi_x,phi_y]=gradient(phi);
s=sqrt(phi_x.^2 + phi_y.^2);
a=(s>=0) & (s<=1);
b=(s>1);
ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  % compute first order derivative of the double-well potential p2 in eqaution (16)
dps = ps./(s+eps);
[dps_x, dps_y] = gradient(dps);
f = (dps_x.*phi_x + dps_y.*phi_y) + 4*dps.*del2(phi);

%heaviside function
% function h = Heaviside(x,epsilon)    
% h=0.5*(1+(2/pi)*atan(x./epsilon));
function h = Heaviside(x,epsilon)
h = x>0;

% Make a function satisfy Neumann boundary condition
function g = NeumannBoundCond(f)
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  

%div(\nabla\phi / |\nabda \phi|)
function k = div_norm(u)
% compute curvature for u with central difference scheme
[ux,uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-10);
Nx = ux./(normDu+eps);
Ny = uy./(normDu+eps);
k = divergence(Nx, Ny);

%dh/dx
function f = Dirac(x, epsilon)    
f=(epsilon/pi)./(epsilon^2+x.^2);

function sparseRecon = sparse_recon(dmDict, sparseCoef)
nSample = size(dmDict,2);
sparseRecon = zeros(size(dmDict{1}));
for i = 1:nSample
    sparseRecon = sparseRecon + dmDict{i}*sparseCoef(i);
end

function u = post_process(u, peakX, peakY)
BW = u>0;
c0 = 5;
[nRow, nCol] = size(u);
CC = bwconncomp(BW,8);
linearPeak = (peakX-1)*nRow+peakY;
nComponent = size(CC.PixelIdxList,2);
biasThreshold = 3;
for iComponent = 1:nComponent
    linearIdx = CC.PixelIdxList{iComponent};
    matX = ceil(linearIdx/nRow);
    matY = linearIdx-(matX-1)*nRow;
%     if isempty(find(linearIdx>=linearPeak-linearThd & linearIdx<=linearPeak+linearThd, 1))
    biasDist = sqrt((matX-peakX).^2+(matY-peakY).^2);
    if isempty(find(biasDist <= biasThreshold, 1))   
        
        u(matY,matX) = -c0;
    end
end

function g = conv2FFTDiv(h, f)

sh = size(h);
sf = size(f);

% zero pad the input signals
fm = zeros(sf+2*(sh-1), class(f));
o = sh-1;
fm( o(1)+(1:size(f,1)), o(2)+(1:size(f,2)) ) = f;

h_zp = zeros(size(fm), class(h));
h_zp(1:size(h,1), 1:size(h,2)) = h;

% compute the convolution in frequency
F = fft2(fm);
H = fft2(h_zp);
Y = F./H;

% back to spatial domain
g = real( ifft2(Y) );

% remove padding
o = floor(1.5*size(h))-1;
g = g( o(1)+(1:size(f,1)), o(2)+(1:size(f,2)) );
return