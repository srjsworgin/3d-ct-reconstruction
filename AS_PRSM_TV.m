%==============================================================================================
% AS_PRSM for 3D TV Image Reconstruction
%==============================================================================================
%==============================================================================================
function [xr,EE, TM]=AS_PRSM_TV(x0,y,TT,num,paras)
    
    x_sum = zeros(size(x0(:)));
    
for T = 1:num
    %============================================================
    xs=paras.xs; % ground truth
    N=size(x0);  % the dimension of 3D image
    n = length(y);
    Amx=paras.Amx; 
    Atmx=paras.Atmx;
    paraA=paras.para;
    paraW.nx=N(1); 
    paraW.ny=N(2);
    paraW.nz=N(3);
    %============================================================
    %%%% set parameter %%%%%%%%%%%%
    maxit = 2000;
    alpha = 6e-2;
    beta =8e-11;
    s = 1.95;
    lambda = 2e-2;%1e-2/n;   %lambda/sigma =1e-5;
    sigma = 1e-8; %lambda/1e-3;
    rho = 0.20;
    mk = 40;
    beta_t = 0.95;
    gamma_t = 1e-2;
    erg_k = 0;
    %============================================================
    paraA.nv_block = uint32(4);
    Batch = single(length(paraA.y_det)*length(paraA.z_det)*paraA.nv_block);
    md = uint32(668/paraA.nv_block);
    %============================================================
    %%% initialize% %%%%%%%%%
    x = x0(:);
    xtilde = x0(:);
    h = x0(:);
    u = cell(3,1);
    z = cell(3,1);
    Par1 = cell(3,1);
    Par2 = cell(3,1);
    Par3 = cell(3,1);
    Par4 = cell(3,1);

    
    for i = 1:3
        
        u{i} = Wxs(x,paraW,i);
        z{i} = u{i};
        Par1{i} = zeros(N(1),N(2),N(3),'single');
        Par2{i} = zeros(N(1),N(2),N(3),'single');
        Par3{i} = zeros(N(1),N(2),N(3),'single');
        Par4{i} = zeros(N(1),N(2),N(3),'single');
        
    end

    
    EE(1,T) = norm(x-xs(:));
    TM(1,T) = 0;
    count = 2;
    %============================================================     
    tic;     
    %============================================================
    for out_iter = 1 : maxit
        
           
        xold = x;
        %%%% update h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:3
            Par1{k} = z{k} - beta*(Wxs(x,paraW,k) - u{k});
        end
        h = -Wtxs(Par1);

        %%%%% set parameter in inner loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nb=1;
        RD = zeros(1,md*nb);
        
        for j = 1:nb
            idx = randperm(md);
            RD(1,(j - 1)*md + 1:j*md) = idx;
        end
        
        for iner_iter =1:mk
           
            
            %%%% update xhat %%%%%%%%%%%%%%%%%%%%%%%%%
            
            xhat = beta_t * xtilde + (1 - beta_t) *x;
            
            %%%% compute stochastic gradient %%%%%%%%%%
            
            paraA.I1 = RD(1,iner_iter);
            paraA.I2 = RD(1,iner_iter);
            sgd = Atmx(Amx(xhat,paraA) - y,paraA);
            
            %%% update xtilde %%%%%%%%%%%%%%%%%%%%%%%%
            temp1 = gamma_t * sigma * xtilde + rho * x - sgd/Batch - h;
            xtilde = temp1./(sigma*gamma_t + rho);
            
            %%% update x %%%%%%%%%%%%%%%%%%%%%%%%
            
            x = beta_t * xtilde + (1 - beta_t)*x;
            
           
            
        end
        
%         max(x)
%         min(x)
        
        
        
%         %%%%%% update rho_k %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         diff_x = x - xold ;
%         deta_1k  = norm( diff_x )^2;  
%         diff_x1 = Wxs(diff_x,paraW,1);
%         diff_x2 = Wxs(diff_x, paraW, 2);
%         diff_x3 = Wxs(diff_x, paraW, 3);
%         deta_2k  =  norm(diff_x1(:))^2 + norm(diff_x2(:))^2 + norm(diff_x3(:))^2;
%         if ( (rho*deta_1k - beta*deta_2k) < 0 )
%             rho_min = eta * rho_min;
%         end
%         rho_BB= beta* deta_2k/deta_1k;
%         rho   = max(rho_min, rho_BB);
         
        %%%% update z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        
        
        for k = 1:3
            
            Par4{k} = Wxs(x,paraW,k) - u{k};
            tempz{k} = z{k} - alpha*beta*(Par4{k});
            Par2{k} = u{k} - (tempz{k} - s*beta*(Par4{k}))/sigma;
            
        end
        
        oldu = u;
        
        %%%% update u %%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        u = Prox(Par2, lambda/sigma);

%         max(max(max(abs(u{1}))))
%         max(max(max(abs(u{2}))))
%         max(max(max(abs(u{3}))))
        
       %%%% update z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       for k =1:3
           
           Par3{k} = s* Wxs(x,paraW,k) + (1-s)*oldu{k} - u{k};
           z{k} = tempz{k} - beta * Par3{k};
           
       end
   
      
     %========================================================================================================
      % Update of error and time.
     %========================================================================================================
        tim1  = toc;
        
        if  tim1 > (TT/1.002) 
            erg_k = erg_k + 1 ;
            alpha_k = 1/erg_k;       
            x_ag = (1- alpha_k)*x_ag + alpha_k*x ;    
        else
            x_ag = x;
      
       end  


       
        ee = norm(x_ag - xs(:));
        fprintf('time:%f,value:%f\n',tim1,ee);
        EE(count,T) = ee;          %norm(x_ag - xs(:));
        TM(count,T) = tim1;
        %PSNR(count, 1) = pps;
        count = count + 1;
        
        if tim1 >= TT
            break;
        end

        

       
        
        %========================================================================================================
    end
    
    x_sum = x_sum + x_ag;
    
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


function ht=HardT(ht,c)
for i = 1:3
ht{i} = ht{i} .* (ht{i} > sqrt(2*c));
end
end

function s = Prox(s,c)
for i = 1:3
    s{i} = sign(s{i}).*(max(abs(s{i}) - c,0));
end
end





