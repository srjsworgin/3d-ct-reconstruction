%==============================================================================================
%==============================================================================================
function [xr,EE,TM]=SARAH_ADMM(x_0,y,TT,num,paras)

    %==============================================================================================
 x_sum = zeros(size(x_0(:)));
 
for T = 1:num

    %============================================================
    xs=paras.xs;
    N=size(x_0);
    n = length(y);
    Amx=paras.Amx;
    Atmx=paras.Atmx;
    paraA=paras.para;
    paraW.nx=N(1);
    paraW.ny=N(2);
    paraW.nz=N(3);
    
    %============================================================
    maxoutit = 1e4;
    max_innerit =2;

    sigma = 1;
    tau = 60;
    beta = 5e-11;
    lambda = 1/10;
    
    %============================================================
    paraA.nv_block = uint32(4);
    Batch = single(length(paraA.y_det)*length(paraA.z_det)*paraA.nv_block);
    md = uint32(668/paraA.nv_block);
    
    %============================================================
    
    x_s = x_0(:);
    count = 2;
    EE(1,T) = norm(x_s-xs(:));
    TM(1,T) = 0;
    
    u = cell(3,1);
    z = cell(3,1);
    Par1 = cell(3,1);
    Par2 = cell(3,1);
    Par3 = cell(3,1);
    Par4 = cell(3,1);
    
    for i = 1:3
        
        z{i} = zeros(N(1),N(2),N(3),'single');
        u{i} = Wxs(x_s,paraW,i);
        Par1{i} = zeros(N(1),N(2),N(3),'single');
        Par2{i} = zeros(N(1),N(2),N(3),'single');
        Par3{i} = zeros(N(1),N(2),N(3),'single');
        Par4{i} = zeros(N(1),N(2),N(3),'single');
        
    end
    
    tic;
    
 for s = 1:maxoutit
     
     x0 = x_s;
    
     paraA.I1 = 1;
     paraA.I2 = md;
     v = Atmx(Amx(x0,paraA) - y,paraA)/n; 
     

    
    %%%%%%%%%%%%%% use gd_full update%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for i = 1:3
         Par1{i} = Wxs(x0,paraW,i) + u{i};
     end
     z = HardT(Par1,lambda/beta);
     
    
    
     for i = 1:3
         Par2{i} = Wxs(x0,paraW,i) - z{i}+ u{i};
     end
     x = x0 - tau*(v + beta.* Wtxs(Par2));
     
  
    
     for i = 1:3
         u{i} = u{i} +  sigma.*(Wxs(x,paraW,i) - z{i});
     end
   
    %============================================================   
       
     %md = floor(allag/bsize); %number of batch
     nb=1;
     RD = zeros(1,md*nb);
   

        for j = 1:nb
            idx = randperm(md);
            RD(1,(j - 1)*md + 1:j*md) = idx;
        end

    
   
    for it = 1 : max_innerit
        
        xold = x0;

        
       %%%% update z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       for k = 1:3
           Par3{k} = Wxs(x,paraW,k) + u{k};
       end
       z = Prox(Par3,lambda/beta);
       
%         max(max(max(abs(z{1}))))
%         max(max(max(abs(z{2}))))
%         max(max(max(abs(z{3}))))
 
     
             x0 = x;
       %%%%% compute sto grdient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
                       
            paraA.I1 = RD(1,it);
            paraA.I2 = RD(1,it);
            sgd1 = Atmx(Amx(x,paraA) - y,paraA);
            sgd2 = Atmx(Amx(xold,paraA) - y,paraA);
            v = (sgd1 - sgd2)/Batch + v;
         
            
       %%%%%%%%%%%%%% update x %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            for k = 1:3
                Par4{k} = Wxs(x,paraW,k) - z{k}+ u{k};
            end
            x = x - tau*(v + beta.* Wtxs(Par4));
            
%             max(x)
%             min(x)
       %%%%%%%%%%%%%% update u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1:3
                u{k} = u{k} +  sigma.*(Wxs(x,paraW,k) - z{k});
            end
            
    end 
    
    tim1 = toc;
    x_s = x;
   
   
    ee = norm(x - xs(:));
    fprintf('time:%f,value:%f\n',tim1,ee);
    %pps = 10*log10(256^2*64./(ee.^2))
    
    EE(count,T) = ee;
    TM(count, T) = tim1;
   % PSNR(count, 1) = pps;
    
       if tim1 >= TT
           break;
       end
       

       count = count + 1;

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
       












