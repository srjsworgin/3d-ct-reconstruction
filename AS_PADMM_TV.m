%==============================================================================================
% SG-GAFBA for 3D TV Image Reconstruction
%==============================================================================================
function [xr,EE, TM]=AS_PADMM_TV(x0,b,TT,num,paras)
    
    x_sum = zeros(size(x0(:)));
    
for T = 1:num
    %============================================================
    xs=paras.xs; % ground truth
    N=size(x0);  % the dimension of 3D image 256*256*64
    Amx=paras.Amx; 
    Atmx=paras.Atmx;
    paraA=paras.para;
    paraW.nx=N(1);  % 256
    paraW.ny=N(2);  % 256
    paraW.nz=N(3);  % 64
    %============================================================
    %%%% set parameter %%%%%%%%%%%%
    maxit = 1e4;
%   alpha = 0.5;   
%   mu = 0.0;    
%   bound = (alpha-(mu-1-mu^2)*(1-alpha)^2+...
%          sqrt(((mu-1-mu^2)*(1-alpha)^2+alpha)^2+4*alpha*(1-alpha)^2))/2;
%   Lip = 1e-2;
    sigma = 4e5;
    beta = 1e-6;
    gamma = 1;
    tau = 0.1;
    pho = 1e-5;%0.02
    pho_min=1e-5;
    lambda = 1e-1; 
    mk = 10;         % inner loop     
    etak = 1e7;      % theta_k, inner-loop
    const_eta_k = 2/etak;  % 2/theta_k,inner-loop
%     beta_t=0.95;
%     gamma_t=1e-2;
    %============================================================
    paraA.nv_block = uint32(4);
    Batch = single(length(paraA.y_det)*length(paraA.z_det)*paraA.nv_block);
    md = uint32(668/paraA.nv_block);  % the number of viewers is 668
    %============================================================
    %%% initialize% %%%%%%%%%
    x = x0(:);
    xtilde = x0(:);
    h = x0(:);
    y = cell(3,1);  % y
    z = cell(3,1);
    Par1 = cell(3,1);
    Par2 = cell(3,1);
    Par3 = cell(3,1);
    Par4 = cell(3,1);
    
    for i = 1:3       
        y{i} = Wxs(x,paraW,i);
        z{i} = y{i};
        Par1{i} = zeros(N(1),N(2),N(3),'single'); 
        Par2{i} = zeros(N(1),N(2),N(3),'single');
        Par3{i} = zeros(N(1),N(2),N(3),'single');
        Par4{i} = zeros(N(1),N(2),N(3),'single');
    end   
    EE(1,T) = norm(x-xs(:));
    TM(1,T) = 0;
    count = 2;
    %============================================================     
    tic = cputime;     
    %============================================================
    for out_iter = 1 : maxit   
        %%%% update h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:3
            Par1{k} = z{k} - beta*(Wxs(x,paraW,k) - y{k});
        end
        h = -Wtxs(Par1);

        %%%%% set parameter in inner loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nb=1;
        RD = zeros(1,md*nb);   % 1, 167; md = 167
        
        for j = 1:nb
            idx = randperm(md);
            RD(1,(j - 1)*md + 1:j*md) = idx;
        end        
        oldx = x;  
        for iner_iter =1:mk
             beta_t   = 2/(iner_iter+1); 
             gamma_t  = const_eta_k/iner_iter;             
            %%%% update xhat %%%%%%%%%%%%%%%%%%%%%%%%%
            xhat = beta_t * xtilde + (1 - beta_t) *x;          
            
            %%%% compute stochastic gradient %%%%%%%%%%           
            paraA.I1 = RD(1,iner_iter);
            paraA.I2 = RD(1,iner_iter);
            sgd = Atmx(Amx(xhat,paraA) - b,paraA);
            
            %%% update xtilde %%%%%%%%%%%%%%%%%%%%%%%%
            temp1 = gamma_t * sigma * xtilde + pho*oldx - sgd/Batch - h;   %%%%%%%%  y? %%为什么是-，为什么*tau，h代表啥
            xtilde = temp1./(gamma_t * sigma + pho);
            
            %%% update x %%%%%%%%%%%%%%%%%%%%%%%%            
            x = beta_t * xtilde + (1 - beta_t)*x;
        end
        
        %         %%%%%% update rho_k %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         diff_x = x - oldx ;
%         deta_1k  = norm( diff_x )^2;  
%         diff_x1 = Wxs(diff_x,paraW,1);
%         diff_x2 = Wxs(diff_x, paraW, 2);
%         diff_x3 = Wxs(diff_x, paraW, 3);
%         deta_2k  =  norm(diff_x1(:))^2 + norm(diff_x2(:))^2 + norm(diff_x3(:))^2;
%         if ( (pho*deta_1k - beta*deta_2k) < 0 )
%             pho_min = eta * pho_min;
%         end
%         pho_BB= beta* deta_2k/deta_1k;
%         pho   = max(pho_min, pho_BB);

        %%%% update y %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:3        
            Par2{k} = Wxs(x,paraW,k)+tau*y{k} - z{k}/beta;  
        end   
        y = Prox(Par2, lambda/(1+tau));
        
       %%%% update z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
       for k =1:3   
           Par3{k} = Wxs(x,paraW,k) - y{k};
           z{k} = z{k} - gamma * beta* Par3{k};
       end
            
     %% Update of error and time.
        tim1 = cputime - tic;   
        ee = norm(x - xs(:));
        fprintf('time:%f,value:%f\n',tim1,ee);
        EE(count,T) = ee;          %norm(x_ag - xs(:));
        TM(count,T) = tim1;
        count = count + 1;        
        if tim1 >= TT
            break;
        end      
        %========================================================================================================
    end    
    x_sum = x_sum + x; 
end
xr = x_sum/num;
end


function res = Wxs(image,paraW,option)
image=reshape(image,[paraW.nx,paraW.ny,paraW.nz]);
switch(option)
    case 1
        res = image([2:end,end],:,:) - image;
    case 2
        res = image(:,[2:end,end],:) - image;
    case 3
        res = image(:,:,[2:end,end]) - image;
end
end


function res = Wtxs(y)
res = zeros(size(y(1)));

yy1 = y{1};
par1 = yy1([1,1:end-1],:,:) - yy1;
par1(1,:,:) = -yy1(1,:,:);
par1(end - 1,:,:) = yy1(end - 1,:,:);

yy2 = y{2};
par2 = yy2(:,[1,1:end-1],:) - yy2;
par2(:,1,:) = -yy2(:,1,:);
par2(:,end - 1,:) = yy2(:,end - 1,:);

yy3 = y{3};
par3 = yy3(:,:,[1,1:end-1]) - yy3;
par3(:,:,1) = -yy3(:,:,1);
par3(:,:,end - 1) = yy3(:,:,end - 1);

res = par1 + par2 + par3;
res = res(:);
end

function s = Prox(s,c)
for i = 1:3
    s{i} = sign(s{i}).*(max(abs(s{i}) - c,0));
end
end


