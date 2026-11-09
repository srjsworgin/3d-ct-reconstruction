# 3d-ct-reconstruction

三维 CT 重建实验代码，使用 MATLAB、C++ 和 CUDA。压缩包中的源码、主实验数据和编译文件已按原目录结构展开在本仓库中，可直接点击浏览。

## 下载

- **[下载仓库文件 ZIP](https://github.com/srjsworgin/3d-ct-reconstruction/archive/refs/heads/main.zip)**：包含源码、Prostate.mat、MEX 和随附编译依赖，不包含超大文件 y.mat。也可点击绿色 Code → Download ZIP。
- **[单独下载 y.mat（约 280 MB）](https://github.com/srjsworgin/3d-ct-reconstruction/releases/download/v1.0.0/y.mat)**：下载后放到 `X_ray_transform_3d/y.mat`，恢复完整原目录。
- [下载原始完整 RAR 备份](https://github.com/srjsworgin/3d-ct-reconstruction/releases/download/v1.0.0/3D-CT.reconstruction.problem.code.rar)。

## 文件说明

- `main_3DTV_test.m`：主实验入口。
- `AS_PADMM_TV.m`、`AS_PRSM_TV.m`、`PDHG.m`、`SARAH_ADMM.m`、`stoADMM.m`：算法实现。
- `fig_prsm.m`：绘图脚本。
- `Prostate.mat`：主实验使用的数据。
- `X_ray_transform_3d/`：投影算子源码、MEX 文件和编译脚本。

`y.mat` 超过 GitHub 普通 Git 仓库的单文件限制，因此单独放在 Release 中；上述链接可直接下载，无需进入 Release 页面。主入口会生成投影数据；子目录中的 `test_3d.m` 则引用 y.mat，并含原有 phantom_3d.mat 引用，运行前需检查所需输入。

## 运行环境

在 MATLAB 中进入仓库目录，按本机环境检查 CUDA 和 MEX 依赖。原编译脚本含 CUDA 9.0 和 Visual Studio 14.0 路径，需要按本机环境调整。

原始代码未修改，本次上传未执行实验或验证数值结果。第三方作者署名和版权说明保留在原文件中。
