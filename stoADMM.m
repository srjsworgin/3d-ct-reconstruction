%==============================================================================================
%==============================================================================================
function [xr,EE, TM]=stoADMM(x0,y,TT,num,paras)

     x_sum = zeros(size(x0(:)));
     
for T=1:num
    %============================================================
    xs=paras.xs;
    N=size(x0);
    n = length(y);
    Amx=paras.Amx;
    Atmx=paras.Atmx;
    paraA=paras.para;
    paraW.nx=N(1);
    paraW.ny=N(2);
    paraW.nz=N(3);
    %============================================================
    maxit = 1e4;
    sigma = 0.95;
    tau = 60;
    lambda = 1/20;
    beta = 1e-11;
    
    %============================================================
    paraA.nv_block = uint32(30);
    Batch = single(length(paraA.y_det)*length(paraA.z_det)*paraA.nv_block);
    md = uint32(668/paraA.nv_block);
    %============================================================
    x = x0(:);
    u = cell(3,1);
    z = cell(3,1);
    Par1 = cell(3,1);
    Par2 = cell(3,1);

    
    for i = 1:3
        
        z{i} = zeros(N(1),N(2),N(3),'single');
        u{i} = Wxs(x,paraW,i);
        Par1{i} = zeros(N(1),N(2),N(3),'single');
        Par2{i} = zeros(N(1),N(2),N(3),'single');
        
    end

    
    EE(1,T) = norm(x-xs(:));
    TM(1,T) = 0;
    count = 2;
    %============================================================     
    tic;
    nb=1;
    RD = zeros(1,md*nb);
   
 
        for j = 1:nb
            idx = randperm(md);
            RD(1,(j - 1)*md + 1:j*md) = idx;
        end
     
    %============================================================
    for i = 1 : maxit

        
        
        %%%% update z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:3
            Par1{k} = Wxs(x,paraW,k) + u{k};
        end
       z = Prox(Par1,lambda/beta);
       
%         max(max(max(abs(z{1}))))
%         max(max(max(abs(z{2}))))
%         max(max(max(abs(z{3}))))

        
        %%% compute stochastic gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
            if mod(i, md) == 0
                bdx = md;
            else
                bdx = mod(i,md);
            end
            
            paraA.I1 = RD(1,bdx);
            paraA.I2 = RD(1,bdx);
            sgd = Atmx(Amx(x,paraA) - y,paraA);
        
              
     %%%%%%%%%%%%%% update x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:3
                Par2{k} = Wxs(x,paraW,k) - z{k}+ u{k};
            end
        x = x - tau*(sgd/Batch + beta.* Wtxs(Par2));
        
%         max(x)
%         min(x)
     
     %%%%%%%%%%%%%% update u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     for k = 1:3
         u{k} = u{k} +  sigma.*(Wxs(x,paraW,k) - z{k});
     end
      
     %========================================================================================================
      % Update of error and time.
     %========================================================================================================
        tim1  = toc;
        ee = norm(x - xs(:));
        fprintf('time:%f,value:%f\n',tim1,ee);
        %pps = 10*log10(256^2*64./(ee.^2))
        
        EE(count,T) = ee;
        TM(count, T) = tim1;
        %PSNR(count, 1) = pps;
        
        if tim1 >= TT
           break;
        end

        count = count + 1;
       
       
        
        %========================================================================================================
    end
x_sum = x_sum + x;
end
xr = x_sum/T;
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









