This class WGStoGCJ implements geodetic coordinate transform between geodetic coordinate in WGS84 coordinate system and geodetic coordinate in GCJ-02 coordinate system with high precision.

<font color = "#2A80ff">Supports 4 ways to convert coordinates form GCJ-02 coordinate system to WGS84 coordinate system:</font>
- Simple linear iteration
- Analytical derivation
- Numerical derivation
- Automatic derivation(needed Ceres library)

**Notice:**
1. Use C++ standard 17.
2. Ceres is not a mandatory option for this project, and we suggest to set 'USE_CERES = OFF'. If set check 'USE_CERES', you should set the following librarys path:
[Ceres](https://github.com/ceres-solver/ceres-solver "Ceres")  [Glog](https://github.com/google/glog "Glog") [GFlags](https://github.com/gflags/gflags "GFlags") [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen")

- **Author**: kikkimo
- **Email**: [fywhu@outlook.com](fywhu@outlook.com "fywhu@outlook.com")
- **Details**: https://blog.csdn.net/gudufuyun/article/details/106738942
- School of Remote Sensing  and Information Engineering,WuHan University,Wuhan, Hubei, P.R. China. 430079
- Copyright (c) 2020.  All rights reserved.

> This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.


<br><br>
这个类实现了WGS84大地坐标与GCJ-02大地坐标的相互转换，具有较高的转换精度，可通过迭代实现误差1mm以内的逆转换。
<font color = "#2A80ff">支持4种方式将GCJ-02坐标转换为WGS84坐标:</font>
- 简单线性迭代
- 解析求导
- 数值求导
- 自动求导(需借助Ceres库)


**注意事项：**
1、需使用C++17标准
2、Ceres并非本工程的必选项，使用CMake创建工程时，建议不要勾选USE_CERES。如果勾选了USE_CERES，需要将以下库引入至本工程：
[Ceres](https://github.com/ceres-solver/ceres-solver "Ceres")  [Glog](https://github.com/google/glog "Glog") [GFlags](https://github.com/gflags/gflags "GFlags") [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen")

- **作者**：kikkimo
- **邮箱**：[fywhu@qq.com](fywhu@outlook.com "fywhu@qq.com")
- **原理**：https://blog.csdn.net/gudufuyun/article/details/106738942
- **地址**：武汉大学遥感信息工程学院，中国湖北省武汉市，430079
- Copyright (c) 2020.  All rights reserved.
