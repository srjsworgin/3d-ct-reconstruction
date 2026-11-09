
%===========================================================================
%===========================================================================
close all;
clear;
clc;
N = [256,256,64]; 	
load Prostate.mat
addpath([pwd '/X_ray_transform_3d']);
xs = Image;   % 256*256*64
xs =xs(:)/max(xs(:));
xs = single(xs); % 4194304*1

sigma = single(max(xs)*0.01);
xsn = xs + sigma*randn(size(xs),'single');
N_iter = 3*N(1);
para = set_geo_para_3d;
tic;y = Ax3_cone_mf(xs,para);toc; % 131334144*1
y = y + sigma*randn(size(y),'single');
paras.xs = xs;
x0 = zeros(N,'single');
paras.para = para;

% ===========================================================================
paras.Amx=@Ax3_cone_mf; 
paras.Atmx=@Atx3_cone_mf;

%===========================================================================
paras.mu=single(1);
%===========================================================================
params.tol =0.01;
TT=2000; 
tol=1e-6;
num = 1;

disp('Computation of ASP_ADMM') 

[x_ASP, EE_ASP, TM_ASP] = AS_PADMM_TV(x0,y,TT,num,paras);
x_ASP=reshape(x_ASP,N);

disp('Computation of PDHG') 

[x_PDHG, EE_PDHG, TM_PDHG] = PDHG(x0,y,TT,num,paras);
x_PDHG=reshape(x_PDHG,N);

disp('Computation of AS-PRSM') 

[x_AS, EE_AS, TM_AS] = AS_PRSM_TV(x0,y,TT,num,paras);
x_AS=reshape(x_AS,N);

disp('Computation of SARAH_ADMM')
 
[x_SARAH, EE_SARAH, TM_SARAH] = SARAH_ADMM(x0,y,TT,num,paras);
x_SARAH=reshape(x_SARAH,N);

disp('Computation of stoADMM')

[x_sto, EE_sto, TM_sto] = stoADMM(x0,y,TT,num,paras);
x_sto=reshape(x_sto,N);

ASP_PS =  10*log10(256^2*64./(EE_ASP).^2);
AS_PS =  10*log10(256^2*64./(EE_AS).^2);
sto_PS = 10*log10(256^2*64./(EE_sto).^2);
SARAH_PS = 10*log10(256^2*64./(EE_SARAH).^2);
PDHG_PS =  10*log10(256^2*64./(EE_PDHG).^2);


