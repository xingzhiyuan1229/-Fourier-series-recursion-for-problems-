# 数值方法公式与代码对照解读（第一部分）

本文档依据论文 *Approximate Fourier series recursion for problems involving temporal fractional calculus*（论文翻译.md）中 **第2节（Volterra第二类方程）** 与 **第3节（Fourier级数近似）** 的全部公式，逐一找到库中对应的最底层计算代码，以 **代码在前、公式在后** 的格式进行分段解读。

> 代码仓库地址：[https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-)  
> 当前分支：`copilot/add-numerical-method-formula-interpretation`  
> 主要源文件：
> - [`fouvol.py`](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py)  
> - [`fouker.py`](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py)

---

## 第2节 Volterra第二类方程（Section 2: Volterra equations of the second kind）

本节建立了数值求解所依据的数学模型，给出了问题的精确表述及测试用精确解和右端函数。

---

### 2.1 Volterra第二类积分方程——公式 (10)

**对应代码（fouvol.py，参数初始化与主循环结构）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L87)

```python
# fouvol.py 第87-99行：参数默认值定义，直接对应公式 (10) 中各符号
# defaults definitions
beingloud     =  0              # level of verbosity
smoothness    =  1              # interpolant continuity
solver        =  1              # default to product rectangle rule
Nlim          =  4              # two-sided sum limits for Fourier Series
Nt            =  4              # number of time steps over (0,T)
N1            =  Nt             # number of time steps over (0,T1)
varphi0       =  0.2            # 对应公式中的 phi_0
alpha         = -0.4            # 对应公式中的 alpha，要求 alpha in (-1,0)
T1            =  0.5
T             =  math.pi-T1     # 对应公式中的 T，最终时间
Tx            =  T+T1           # 扩展时间 Tx = T + T1
tbar          =  0.0            # 对应公式中的 t_bar（t上方的横线）
Nvals         =  10000          # for smooth plots
```

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L154)

```python
# fouvol.py 第154-163行：时间步长和 N1 的计算（离散化参数设置）
# 这些参数直接对应公式 (10) 的离散化框架
dt = T/Nt
N1=int(math.ceil(T1/dt))
# crude hack to get into bigrun.sh tables
T1 = N1*dt
Tx = T + T1
if beingloud > 19:
  print "\nSanity adjustments give ", 
  print('T = %f; Nt = %d; dt = %f; T1 = %f; N1 = %d; Tx = %f' % (T,Nt,dt,T1,N1,Tx) )
```

**论文公式 (10)**

公式 (10) 是本文所求解问题的核心方程——带幂律核的 Volterra 第二类积分方程：寻找 \( u : I := [0, T] \to \mathbb{R} \)，使得

\[
u(t) + \int_{0}^{t} \phi(t - s)\, u(s)\, ds = f(t), \quad t \in [0, T]
\]

其中核函数取幂律形式：

\[
\phi(t) = \phi_0 \cdot (t + \bar{t})^{\alpha}, \quad \alpha \in (-1,\, 0),\quad \bar{t} \geq 0,\quad \phi_0 > 0
\]

代码中：
- \( \phi_0 \) 对应 `varphi0`（默认值 0.2）
- \( \alpha \) 对应 `alpha`（默认值 -0.4，需在 \((-1,0)\) 内）
- \( \bar{t} \) 对应 `tbar`（默认值 0.0，即无正则化偏移）
- \( T \) 对应 `T`，\( \Delta t \) 对应 `dt = T/Nt`，\( N_t \) 对应 `Nt`

当 \( \bar{t} = 0 \) 时，核退化为标准弱奇异幂律核，与分数阶微积分直接关联。

---

### 2.2 精确解的设定

**对应代码（fouvol.py，第29-30行：精确解函数）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L29)

```python
# fouvol.py 第29-30行：精确解 u(t)
# exact solution, u(t) - assumes that tbar = 0
def ut(t,alpha):
  return t**(-alpha)
```

**精确解设定**

为了在数值实验中观测收敛阶，论文取如下精确解（假设 \( \bar{t} = 0 \)）：

\[
u(t) = t^{-\alpha}
\]

代码中函数 `ut(t, alpha)` 直接实现此公式，返回 \( t^{-\alpha} \)（即 `t**(-alpha)`）。由于 \( \alpha \in (-1,0) \)，有 \( -\alpha \in (0,1) \)，故 \( u(t) = t^{-\alpha} \) 在 \( t>0 \) 时为有界正函数。

---

### 2.3 卷积积分的解析计算与公式 (11) 的推导基础

**对应代码（fouvol.py，第33-34行：右端函数）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L33)

```python
# fouvol.py 第33-34行：右端函数 f(t)，实现公式 (11)
# RHS function, f(t) - assumes that tbar = 0
def ft(t,varphi0,alpha):
  return t**(-alpha) + varphi0*math.pi*(-alpha)*t/math.sin(-math.pi*alpha)
```

**卷积积分的计算推导**

将精确解 \( u(t) = t^{-\alpha} \) 代入公式 (10) 左端，需要计算卷积积分：

\[
\int_{0}^{t} \phi(t-s)\, u(s)\, ds = \phi_0 \int_{0}^{t} (t-s)^{\alpha}\, s^{-\alpha}\, ds
\]

令 \( s = tu \)，则上式变为：

\[
= \phi_0 \int_{0}^{1} (t - tu)^{\alpha}\, (tu)^{-\alpha}\, t\, du = \phi_0\, t \int_{0}^{1} (1-u)^{\alpha}\, u^{-\alpha}\, du
\]

括号内的积分恰好是 Beta 函数：

\[
\int_{0}^{1} (1-u)^{\alpha}\, u^{-\alpha}\, du = B(1+\alpha,\, 1-\alpha) = \Gamma(1+\alpha)\,\Gamma(1-\alpha)
\]

利用 Gamma 函数的反射公式 \( \Gamma(z)\,\Gamma(1-z) = \frac{\pi}{\sin(\pi z)} \)，取 \( z = 1+\alpha \)（其中 \( \Gamma(1+\alpha) = \alpha\,\Gamma(\alpha) \)）：

\[
\Gamma(1+\alpha)\,\Gamma(1-\alpha) = \alpha\,\Gamma(\alpha)\,\Gamma(1-\alpha) = \frac{\alpha\,\pi}{\sin(\alpha\pi)}
\]

故卷积积分等于：

\[
\phi_0 \int_{0}^{t} (t-s)^{\alpha}\, s^{-\alpha}\, ds = \phi_0 \cdot t \cdot \frac{\alpha\,\pi}{\sin(\alpha\pi)}
\]

**论文公式 (11)**

将上述卷积结果代入 \( f(t) = u(t) + (\phi * u)(t) \)，得到右端函数：

\[
f(t) = t^{-\alpha} + \phi_0\, t^{\alpha} \frac{\pi}{\sin(\alpha\pi)}
\]

> **注意**：论文中写为 \( \phi_0 t^{\alpha}\pi/\sin(\alpha\pi) \)，但依据 Beta 函数严格推导，正确结果应为 \( \phi_0 \cdot t \cdot \frac{\alpha\pi}{\sin(\alpha\pi)} \)，代码实现的正是后者。代码返回值为：
>
> \[
> t^{-\alpha} + \frac{\phi_0 \cdot \pi \cdot (-\alpha) \cdot t}{\sin(-\pi\alpha)} = t^{-\alpha} + \frac{\phi_0 \cdot \alpha\pi\, t}{\sin(\alpha\pi)}
> \]

代码中：
- `t**(-alpha)` 计算 \( t^{-\alpha} \)（精确解）
- `varphi0*math.pi*(-alpha)*t/math.sin(-math.pi*alpha)` 计算卷积项 \( \phi_0 \cdot \alpha\pi t / \sin(\alpha\pi) \)

由于 \( \alpha \in (-1,0) \)，`-alpha > 0`，`math.sin(-math.pi*alpha) = -sin(\pi\alpha) < 0`，故整个第二项为正值，与物理预期一致。

---

### 2.4 数值精度验证代码

**对应代码（fouvol.py，第410-411行：计算误差）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L410)

```python
# fouvol.py 第410-411行：计算数值解与精确解之间的误差
# calc max error
error = abs(ut(T,alpha)-uvals[Nt])
```

**误差量的含义**

该代码计算的是最终时刻 \( t = T \) 处的绝对误差：

\[
\text{error} = \left| u(T) - U_{N_t} \right| = \left| T^{-\alpha} - U_{N_t} \right|
\]

其中 \( U_{N_t} \) 是数值解在最后时间步的值，\( u(T) = T^{-\alpha} \) 是精确解。

---

## 第3节 Fourier级数近似（Section 3: The Fourier series approximation）

本节介绍将 Volterra 积分的"远历史"部分用 Fourier 级数代理的核心思想，包括 Hermite 样条的构造与 Fourier 系数的计算。

---

### 3.1 Volterra积分的分段分割

**对应代码（fouvol.py，第344-351行：solver==3，矩形-Fourier法中的历史分割）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L344)

```python
# fouvol.py 第343-351行：历史积分的分割与 Fourier 贡献
# quadrature over j={i-N1,...,i-2,i-1}
jlo = int(max(1,i+1-N1))
hist = varphi0/(alpha+1)*sum([(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*uvals[j] for j in range(jlo,i)])
# Fourier contribution for far history
if i > N1:
    hist = hist + phistarR[0+Nlim]
    hist = hist + 2*sum([phistarR[j+Nlim] for j in range(-Nlim,0)])
```

**论文中的积分分割公式**

算法的核心思想是将 Volterra 积分分成两段：在 \( [0, T_1] \) 上用标准数值积分处理奇异核，在 \( [T_1, t] \) 上用 Fourier 级数代理。即对 \( t \in [0, T] \)：

\[
\int_{0}^{t} \phi(t-s)\, u(s)\, ds = \int_{0}^{\min\{T_1,\, t\}} \phi(t-s)\, u(s)\, ds + \int_{\min\{T_1,\, t\}}^{t} \phi(t-s)\, u(s)\, ds
\]

在远历史段 \( [0, \min\{T_1, t\}] \) 中用 \( \psi \) 替换 \( \phi \)（当 \( t > T_1 \) 时）：

\[
= \int_{0}^{\min\{T_1,\, t\}} \phi(t-s)\, u(s)\, ds + \int_{\min\{T_1,\, t\}}^{t} \psi(t-s)\, u(s)\, ds
\]

其中 \( \psi \in C^m_{2T_x\text{-ev-per}}(\mathbb{R}) \) 是偶函数且 \( \psi\big|_{[T_1, T]} = \phi\big|_{[T_1, T]} \)。

代码对应关系：
- `hist` 的第一部分（`jlo` 到 `i` 的求和）对应近历史段的数值积分，范围为 \( [t_{i-N_1}, t_i] \)
- `phistarR[0+Nlim]` 和 `phistarR[j+Nlim]` 对应远历史段的 Fourier 级数贡献

---

### 3.2 Hermite样条构造——左端插值（pL 的构造）

**对应代码（fouker.py，第282-310行：get_FourierCoefficientsNew 中左端 Hermite 计算）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py#L282)

```python
# fouker.py 第282-310行：构造左端 Hermite 插值（对应 p_L 在 [0, T1] 上的构造）
def get_FourierCoefficientsNew(tbar,T1,T,Tx,varphi0,alpha,smoothness,Nlim,DDuse,loudness):
  lval, rval = get_HermiteInterceptValues(0.9,0.5,varphi0,T1,T,Tx,tbar,alpha)
  #calculate the Hermite at the left end
  if DDuse == 1:
    # get the function values
    fa = lval; fb = varphi0*(T1+tbar)**alpha;
    # get the derivative values (all zero on left), d[1] = first deriv etc
    da = np.zeros(1+smoothness); db = np.zeros(1+smoothness)
    deriv = alpha*fb/(T1+tbar) # init to first deriv at T1
    for n in range(1,1+smoothness):
      db[n] = deriv
      deriv = (alpha-n)*deriv/(T1+tbar)
    DDmatL, zL = get_DDHermite(0,T1,smoothness,fa,fb,da,db)
  else:    
    a=0.; b=1.
    fL=np.zeros((2*smoothness+2,1)); 
    fL[0] = lval
    fL[1] = varphi0*(T1+tbar)**alpha
    for k in range(2,2*smoothness+2,2):
      fL[k]   = 0.0
      fL[k+1] = fL[k-1] * (alpha-(k-2.0)/2.0) * (T1+tbar)**(-1)
    pL = get_basicHermite(smoothness,fL,T1-0.,loudness)
```

**左端 Hermite 多项式 \( p_L \) 的边界条件**

论文规定左端 Hermite 多项式 \( p_L : [0, T_1] \to \mathbb{R} \)（次数为 \( 2m+1 \)）满足：

\[
p_L^{(k)}(0) = 0, \quad k = 1, 2, \ldots, m
\]

\[
p_L^{(0)}(0) = p_L(0) \text{ 在 } \phi(T_1) \text{ 与切线截距之间选取}
\]

\[
p_L^{(k)}(T_1) = \phi^{(k)}(T_1), \quad k = 0, 1, \ldots, m
\]

代码中 `da = np.zeros(1+smoothness)` 实现了 \( p_L^{(k)}(0) = 0 \)（左端所有导数为零）的条件，`db[n] = deriv` 则根据 \( \phi(t) = \phi_0(t+\bar{t})^\alpha \) 的各阶导数递推计算右端点 \( T_1 \) 处的值：

\[
\phi^{(n)}(T_1) = \alpha(\alpha-1)\cdots(\alpha-n+1)\,\phi_0\,(T_1+\bar{t})^{\alpha-n}
\]

在代码中以迭代形式实现（`deriv = (alpha-n)*deriv/(T1+tbar)`）。

---

### 3.3 Hermite样条构造——右端插值（pR 的构造）

**对应代码（fouker.py，第312-331行：右端 Hermite 计算）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py#L312)

```python
# fouker.py 第312-331行：构造右端 Hermite 插值（对应 p_R 在 [T, Tx] 上的构造）
  #calculate the Hermite at the right end
  if DDuse == 1:
    # get the function values
    fa = varphi0*(T+tbar)**alpha; fb = rval
    # get the derivative values (all zero on right), d[1] = first deriv etc
    da = np.zeros(1+smoothness); db = np.zeros(1+smoothness)
    deriv = alpha*fb/(T+tbar) # init to first deriv at T1
    for n in range(1,1+smoothness):
      da[n] = deriv
      deriv = (alpha-n)*deriv/(T+tbar)
    DDmatR, zR = get_DDHermite(T,T+Tx,smoothness,fa,fb,da,db)
  else:    
    a=0.; b=1.
    fR=np.zeros((2*smoothness+2,1)); 
    fR[0] =     varphi0*(T+tbar)**alpha
    fR[1] = rval
    for k in range(2,2*smoothness+2,2):
      fR[k]   = fR[k-2] * (alpha-(k-2.0)/2.0) * (T+tbar)**(-1)
      fR[k+1] = 0.0
    pR = get_basicHermite(smoothness,fR,Tx-T,loudness)
```

**右端 Hermite 多项式 \( p_R \) 的边界条件**

论文规定右端 Hermite 多项式 \( p_R : [T, T_x] \to \mathbb{R} \)（次数为 \( 2m+1 \)）满足：

\[
p_R^{(k)}(T) = \phi^{(k)}(T), \quad k = 0, 1, \ldots, m
\]

\[
p_R^{(k)}(T_x) = 0, \quad k = 1, 2, \ldots, m
\]

\[
p_R(T_x) \text{ 在 } \phi(T) \text{ 与切线截距之间选取}
\]

代码中 `da[n] = deriv` 按 \( \phi^{(n)}(T) = \alpha(\alpha-1)\cdots(\alpha-n+1)\phi_0(T+\bar{t})^{\alpha-n} \) 递推计算左端（即 \( T \) 处）的导数值，`db = np.zeros(1+smoothness)` 实现 \( p_R^{(k)}(T_x) = 0 \)（右端导数为零）。

---

### 3.4 分差矩阵——Hermite插值的底层实现

**对应代码（fouker.py，第141-166行：get_DDHermite 函数）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py#L141)

```python
# fouker.py 第141-166行：构造 C^m Hermite 插值的分差矩阵
# get C^m Hermite interpolant divided differences matrix with m derivs at a and b
def get_DDHermite(a,b,m,fa,fb,da,db):
  # an empty matrix, and a vector of duplicated x values
  DDmat = np.zeros((2*m+2,2*m+2))
  z     = np.zeros(2*m+2)
  # populate the duplicated values
  for k in range(0,m+1):
    z[k] = a; z[m+1+k] = b;
  # intialize the matrix
  for k in range(0,m+1):
    DDmat[k,k] = fa;
    DDmat[k+m+1,k+m+1] = fb;
  # now fill the matrix up, do the diagonals for duplicated points first
  nfac = 1
  for n in range(1,m+1):        # update super-diagonal zero blocks 
    nfac = n*nfac
    for k in range(0,m+1-n):    # update rows from top and bottom 
      DDmat[k,  k+n]       = da[n]/nfac
      DDmat[m+k+1,m+k+n+1] = db[n]/nfac
  # finish off with the remaining super diagonal sub block, bottom up
  for k in range(m,-1,-1):        # update upper right square block
    for n in range(m+1, 2*m+2):     # update rows 0, 1, ..., 2m+1-n
      DDmat[k,n] = ( DDmat[k+1,n]-DDmat[k,n-1] ) / ( b-a ) # denom is always b-a
  return DDmat, z
```

**Hermite 插值的分差公式**

对于 \( [a, b] \) 上的 \( C^m \) Hermite 插值，利用重复节点 \( z_0 = z_1 = \cdots = z_m = a \) 和 \( z_{m+1} = \cdots = z_{2m+1} = b \) 构造牛顿型插值。重复节点处的高阶分差由导数值给出：

\[
[z_k, z_{k+1}, \ldots, z_{k+n}]f = \frac{f^{(n)}(a)}{n!}, \quad \text{当 } z_k = z_{k+1} = \cdots = z_{k+n} = a
\]

代码中 `DDmat[k, k+n] = da[n]/nfac` 实现了上式（`nfac` 即 \( n! \)），`DDmat[k,n] = (DDmat[k+1,n]-DDmat[k,n-1])/(b-a)` 是标准分差递推公式：

\[
[z_k, z_{k+1}, \ldots, z_{k+n}]f = \frac{[z_{k+1}, \ldots, z_{k+n}]f - [z_k, \ldots, z_{k+n-1}]f}{z_{k+n} - z_k}
\]

---

### 3.5 Hermite插值的求值

**对应代码（fouker.py，第169-175行：get_DDHermiteValue 函数）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py#L169)

```python
# fouker.py 第169-175行：用分差矩阵求 Hermite 插值在点 t 处的值
# get the value of the Divided difference version of the Hermite at t
def get_DDHermiteValue(t,DDmat,m,z):
  total = 0.0
  coeff = 1.0
  for k in range(0,2*m+2):
    total = total+DDmat[0,k]*coeff
    coeff = coeff*(t-z[k])
  return total
```

**牛顿插值公式的求值**

给定分差矩阵，牛顿型插值多项式 \( p(t) \) 的求值公式为：

\[
p(t) = \sum_{k=0}^{2m+1} [z_0, z_1, \ldots, z_k]f \cdot \prod_{j=0}^{k-1}(t - z_j)
\]

代码中 `total` 累加各项，`coeff` 递推维护 \( \prod_{j=0}^{k-1}(t - z_j) \)（初始为 1，每步乘以 `(t - z[k])`）。`DDmat[0, k]` 是第 \( k \) 阶分差（即 \( [z_0, \ldots, z_k]f \)）。

---

### 3.6 复指数形式的Fourier级数展开

**对应代码（fouker.py，第333-352行：Fourier系数的数值计算）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py#L333)

```python
# fouker.py 第333-352行：Fourier 系数的核心计算
# calculate Fourier series
cnR= np.zeros((2*Nlim+1,1)) 
# basic assumption: Im part = 0, Re part symmetric about n=0. Save time.
for n in range(0,Nlim+1):
  if DDuse == 1:
    realpart = quad(lambda t: get_DDHermiteValue(t,DDmatL,smoothness,zL) * math.cos(math.pi*n*t/Tx), 0, T1, limit=80000)
  else:
    realpart = quad(lambda t: get_MappedHornerPolyValue(pL,0,T1,t,2*smoothness+1) * math.cos(math.pi*n*t/Tx), 0, T1, limit=80000)
  cnR[n+Nlim] = 1.0/Tx*realpart[0]
  realpart = quad(lambda t: varphi0*(t+tbar)**alpha * math.cos(math.pi*n*t/Tx), T1, T, limit=80000)
  cnR[n+Nlim] = cnR[n+Nlim] + 1.0/Tx*realpart[0]
  if DDuse == 1:
    realpart = quad(lambda t: get_DDHermiteValue(t,DDmatR,smoothness,zR) * math.cos(math.pi*n*t/Tx), T, Tx, limit=80000)
  else:
    realpart = quad(lambda t: get_MappedHornerPolyValue(pR,T,Tx,t,2*smoothness+1) * math.cos(math.pi*n*t/Tx), T, Tx, limit=80000)
  cnR[n+Nlim] = cnR[n+Nlim] + 1.0/Tx*realpart[0]

for n in range(-Nlim,0):
  cnR[n+Nlim] = cnR[-n+Nlim]
```

**论文中的复指数Fourier级数公式**

论文首先给出 \( \psi \) 的完整复指数 Fourier 展开（无穷级数形式）：

\[
\psi(t) = \sum_{k=-\infty}^{\infty} c_k\, e^{i\pi k t / T_x}
\]

其中 Fourier 系数定义为：

\[
c_k = \frac{1}{2T_x} \int_{-T_x}^{T_x} \psi(t)\, e^{-i\pi k t / T_x}\, dt
\]

代码中 `cnR` 数组存储两侧系数 \( c_{-N}, c_{-N+1}, \ldots, c_{-1}, c_0, c_1, \ldots, c_N \)（索引偏移量为 `Nlim`）。最后的 `for n in range(-Nlim,0): cnR[n+Nlim] = cnR[-n+Nlim]` 利用 \( c_{-k} = c_k \) 的对称性（\( \psi \) 为偶实值函数）进行填充。

---

### 3.7 实余弦形式的Fourier系数计算公式

**对应代码（fouker.py，第336-350行：余弦积分计算 cnR）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py#L336)

```python
# fouker.py 第336-350行：用余弦积分计算 Fourier 系数（利用 psi 的偶函数性质）
for n in range(0,Nlim+1):
  # 在 [0, T1] 段用左端 Hermite 插值 p_L 计算
  realpart = quad(lambda t: get_DDHermiteValue(t,DDmatL,smoothness,zL) * math.cos(math.pi*n*t/Tx), 0, T1, limit=80000)
  cnR[n+Nlim] = 1.0/Tx*realpart[0]
  # 在 [T1, T] 段用真实核 phi(t) = varphi0*(t+tbar)^alpha 计算
  realpart = quad(lambda t: varphi0*(t+tbar)**alpha * math.cos(math.pi*n*t/Tx), T1, T, limit=80000)
  cnR[n+Nlim] = cnR[n+Nlim] + 1.0/Tx*realpart[0]
  # 在 [T, Tx] 段用右端 Hermite 插值 p_R 计算
  realpart = quad(lambda t: get_DDHermiteValue(t,DDmatR,smoothness,zR) * math.cos(math.pi*n*t/Tx), T, Tx, limit=80000)
  cnR[n+Nlim] = cnR[n+Nlim] + 1.0/Tx*realpart[0]
```

**论文中的实余弦Fourier系数公式**

由于 \( \psi \) 是偶实值函数，利用 \( c_{-k} = c_k \in \mathbb{R} \)，可将两侧对称展开化简为单侧余弦级数：

\[
\psi(t) = \sum_{k=0}^{\infty} c_k\, e^{i\pi k t / T_x}
\]

其中 Fourier 系数的实余弦表示为：

\[
c_k = \begin{cases}
\dfrac{2}{T_x} \displaystyle\int_{0}^{T_x} \psi(t)\, \cos\!\left(\dfrac{\pi k t}{T_x}\right) dt, & k \neq 0 \\[10pt]
\dfrac{1}{T_x} \displaystyle\int_{0}^{T_x} \psi(t)\, dt, & k = 0
\end{cases}
\]

> **代码与公式的对应说明**：代码使用系数 `1.0/Tx`（而非 `2.0/Tx`），这是因为代码采用的是**两侧对称系数**（即标准复指数 Fourier 系数）：
>
> \[
> c_k = \frac{1}{2T_x} \int_{-T_x}^{T_x} \psi(t)\, e^{-i\pi kt/T_x}\, dt = \frac{1}{T_x} \int_{0}^{T_x} \psi(t)\, \cos\!\left(\frac{\pi kt}{T_x}\right) dt
> \]
>
> （利用 \( \psi \) 的偶性将 \([-T_x, T_x]\) 上的积分化为 \([0, T_x]\) 上积分的两倍，再除以 \( 2T_x \) 得 \( 1/T_x \) 因子）。
>
> 论文中 \( k \neq 0 \) 时出现的 \( 2/T_x \) 系数是**单侧半范围余弦级数**（one-sided cosine series）的写法，与代码的两侧写法相差正好一个因子 2，但两种表示给出相同的 \( \psi(t) \)，因为两侧求和时 \( k \neq 0 \) 的项各出现两次。

代码将 \( \psi \) 的 \([0, T_x]\) 域分为三段分别数值积分：
- \( [0, T_1] \)：用左端 Hermite 插值 \( p_L \) 替代 \( \psi \)，调用 `get_DDHermiteValue(...DDmatL...)`
- \( [T_1, T] \)：直接用核函数 \( \phi(t) = \phi_0(t+\bar{t})^\alpha \)（代码为 `varphi0*(t+tbar)**alpha`）
- \( [T, T_x] \)：用右端 Hermite 插值 \( p_R \) 替代 \( \psi \)，调用 `get_DDHermiteValue(...DDmatR...)`

数值积分均使用 `scipy.integrate.quad`，精度极高（`limit=80000` 表示允许最多 80000 个子区间的自适应求积分）。

---

### 3.8 Fourier代理截断误差的估算

**对应代码（fouker.py，第208-222行：get_partialerror 函数）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py#L208)

```python
# fouker.py 第208-222行：计算 Fourier 代理的误差（L_inf 范数近似）
def get_partialerror(varphi0,alpha,Nvals,T1,T,tbar,cnR,Nlim,Tx):
  error = 0
  tvals = np.zeros(Nvals+1)
  FSvals = np.zeros(Nvals+1)
  phivals = np.zeros(Nvals+1)
  for i in range(0,Nvals+1,1):
    t        = float(T1+i*(T-T1)/Nvals)
    tvals[i] = t
    FSvals[i] = 0
    # earlier results had range(-Nlim,Nlim,1)
    for k in range(-Nlim,Nlim+1,1):
      FSvals[i] = FSvals[i] + cnR[k+Nlim]*math.cos(math.pi*k*t/Tx)
    phivals[i] = varphi0*(t+tbar)**alpha
    error = max(abs(FSvals[i]-phivals[i]), error)
  return error
```

**Fourier级数截断误差的含义**

该函数估计 \( L_\infty(T_1, T) \) 意义下 Fourier 级数代理对真实核的近似误差：

\[
\text{FSerror} \approx \left\| (I - F_L)\phi \right\|_{L_\infty(T_1, T)} = \max_{t \in [T_1, T]} \left| \phi(t) - \sum_{k=-L}^{L} c_k e^{i\pi kt/T_x} \right|
\]

其中截断级数（代码中用余弦形式，因 \( \psi \) 为偶实值函数）为：

\[
F_L\psi(t) = \sum_{k=-L}^{L} c_k\, \cos\!\left(\frac{\pi k t}{T_x}\right)
\]

代码中 `FSvals[i] = sum(cnR[k+Nlim]*cos(pi*k*t/Tx) for k in range(-Nlim,Nlim+1))` 计算截断 Fourier 级数在 \( t \in [T_1, T] \) 上的值，`phivals[i] = varphi0*(t+tbar)**alpha` 计算真实核值，`error` 追踪二者最大绝对差。

---

## 小结

本部分（第一部分）覆盖了论文第2节和第3节的全部公式，给出了从 Volterra 方程的精确解与右端函数（公式11），到 \( \psi \) 函数的 Fourier 级数表示（复指数形式与余弦形式）的完整代码对照。第二部分（见 `数值方法公式代码解读_Part2.md`）将继续覆盖第4节（矩形-Fourier规则实现）的全部公式，包括乘积矩形规则公式 (12)、递归更新公式 (15) 和最终求解公式 (16)。
