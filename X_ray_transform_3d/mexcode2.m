system('nvcc -c --compiler-options -fPIC Atx3_cone_mf_gpu_new.cu')
system('nvcc -c --compiler-options -fPIC Ax3_cone_mf_gpu_new.cu ')
mex Atx3_cone_mf.cpp  Atx3_cone_mf_gpu_new.o find_l_3d.cpp -lcudart -lcudart -L"/usr/local/cuda-10.0/lib64"
mex Ax3_cone_mf.cpp  Ax3_cone_mf_gpu_new.o -lcudart -L"/usr/local/cuda-10.0/lib64"
