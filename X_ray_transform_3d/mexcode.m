clc;
system('nvcc -c Atx3_cone_mf_gpu_new.cu -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin"')
mex Atx3_cone_mf.cpp  Atx3_cone_mf_gpu_new.obj find_l_3d.cpp -lcudart -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\lib\x64"

system('nvcc -c Ax3_cone_mf_gpu_new.cu -ccbin "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin"')
mex Ax3_cone_mf.cpp  Ax3_cone_mf_gpu_new.obj -lcudart -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v9.0\lib\x64"