## vb2molden.py



```bash

python vb2molden.py <.orb文件> <.xmo文件>

```

(文件名不包含后缀)

该脚本从.orb中读取轨道信息，从.xmo文件中读取分子结构、基组、基函数数量(需要`int=libcint`)

输出与.orb文件名相同的.molden



## molden2gus.py



```bash

python molden2gus.py <.molden文件>

```

借用了pyscf.tools.molden，需要安装pyscf才能使用

该脚本读取.molden文件中的轨道给XMVB提供初猜

输出记录初猜轨道的.gus文件和记录初猜轨道的_gus.molden文件



将Gaussian等软件输出的轨道定域化之后使用此脚本提供初猜



使用方法:

输入`<原子> <轨道>`

例如: `2,4 6-8,9`

将2,4号原子的第6~8,9轨道添加到初猜中(添加4个轨道)



输入`a <轨道>`添加所有原子



输入`r<角度> <轨道1>,<轨道2>`旋转两个轨道(乘对应的2维旋转矩阵)



输入`m<数字> <轨道>`将轨道系数乘以<数字>(可以是负数)



输入`q`退出



可以将输入写到`.txt`文件中使用`< .txt`提供输入



## gus2molden.py



```bash

python gus2molden.py <.xdat文件> <.xmo文件>

```

该脚本读取.xdat文件中的初猜轨道，从.xmo文件中读取分子结构、基组、基函数数量(需要`int=libcint`)

输出与.xdat文件名相同的.molden



用于检查初猜





## sortw.py



```bash

python sortw.py <.xmo文件> <参数(可选)>

```

排序权重/系数

默认参数为w

w: WEIGHTS OF STRUCTURES

l: Lowdin Weights

i: Inverse Weights

r: Renormalized Weights

c: COEFFICIENTS OF STRUCTURES

lc: LOWDIN ORTHOGONALIZED COEFFICIENTS OF STRUCTURES



## no2molden.py



```bash

python no2molden.py <.xmo文件>

```

读取xmvb.no中的自然轨道，从.xmo文件中读取分子结构、基组、基函数数量(需要`int=libcint`)

输出与.xmo文件名相同的.molden



如果怀疑轨道文件不正确可以用Multiwfn的1000 100功能检查波函数是否归一化
