function para=set_geo_para_3d
nx=256;ny=nx;nz=64;
dx=single(250/nx);dy=dx;dz=single(250/nz);
% nv=4;
nv=668;
SO=single(1000);OD=single(500);
na=512;dy_det=single(0.388*1024/na);
nb=384;dz_det=dy_det;
sd_phi=single(2*pi/nv*(0:nv-1));    
sd_z=single(zeros(1,nv)); % for helical scan
y_os=single(0*dy_det);
y_det=single(((-na/2:na/2-1)+0.5)*dy_det+y_os);
z_det=single(((-nb/2:nb/2-1)+0.5)*dz_det);
% normalize so that dx=dy=1 %
scale=dx;
sd_z=sd_z/scale;SO=SO/scale;OD=OD/scale;
y_det=y_det/scale;z_det=z_det/scale;
dy_det=dy_det/scale;dz_det=dz_det/scale;
y_os=y_os/scale;dz=dz/scale;
% parameters (single image frame) %
% note: this version works for multiple images as well by supplying id_v -- which projection views are for each frame%
% load phantom_3d.mat x0 % 3D phantom image
nt=1;
% X0=x0(:);
Id_v=cell(1,nt);
Id_v{1}=uint32(0:nv-1);
id_X=[];Nv=zeros(1,nt);
for i=1:nt
    id_X=[id_X (i-1)*ones(1,numel(Id_v{i}))];
    Nv(i)=numel(Id_v{i});
end
tmp_size=max(Nv);
id_Y=zeros(tmp_size,nt);
for i=1:nt
    id_Y(1:Nv(i),i)=Id_v{i};
end
para=struct('version',[],'GPU',[],...
    'SO',single(SO),'OD',single(OD),'scale',single(scale),'nx',uint32(nx),'ny',uint32(ny),'nz',uint32(nz),'nt',uint32(nt),...
    'sd_phi',sd_phi,'sd_z',sd_z,'y_det',y_det,'z_det',z_det,'cos_phi',single(cos(sd_phi)),'sin_phi',single(sin(sd_phi)),...
    'dy_det',single(dy_det),'y_os',single(y_os/dy_det),'dz_det',single(dz_det),'dz',single(dz),...
    'id_X',uint32(id_X),'id_Y',uint32(id_Y),'Nv',uint32(Nv),'tmp_size',uint32(tmp_size),'nv_block',uint32(4),'I1',uint32(1),'I2',uint32(167));
para.version=uint32(1);para.GPU=uint32(1);