# 数值方法公式与代码对照解读（第二部分）

本文档是 **数值方法公式代码解读_Part1.md** 的续篇，覆盖论文 *Approximate Fourier series recursion for problems involving temporal fractional calculus*（论文翻译.md）中 **第4节（矩形-Fourier规则的实现）** 的全部公式，以 **代码在前、公式在后** 的格式进行分段解读。

> 代码仓库地址：[https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-)  
> 当前分支：`copilot/add-numerical-method-formula-interpretation`  
> 主要源文件：
> - [`fouvol.py`](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py)  
> - [`fouker.py`](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouker.py)

---

## 第4节 矩形-Fourier规则的实现（Section 4: Implementation: the rectangle-Fourier rule）

本节给出将 Fourier 级数代理与乘积矩形规则结合的完整数值格式，包括从基本乘积矩形规则（公式12）到带 Fourier 代理的递归算法（公式15和16）的全部推导和实现。

---

### 4.1 时间步长和 N1 的设置

**对应代码（fouvol.py，第154-163行：离散化参数）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L154)

```python
# fouvol.py 第154-163行：时间步长 dt、N1 和 T1 的精确设置
# sanity checks
dt = T/Nt
N1=int(math.ceil(T1/dt))
# crude hack to get into bigrun.sh tables
#os.system('echo '+N1+' > N1_value.txt')
T1 = N1*dt
Tx = T + T1
if beingloud > 19:
  print "\nSanity adjustments give ", 
  print('T = %f; Nt = %d; dt = %f; T1 = %f; N1 = %d; Tx = %f' % (T,Nt,dt,T1,N1,Tx) )
```

**离散化参数定义**

论文第4节在实现细节中规定（参见论文第302行附近的说明）：

取 \( \Delta t = T / N_t \)，\( t_n = n\Delta t \)，并要求 \( N_1 = T_1 / \Delta t \in \mathbb{N} \)（整数）。由于 \( T_1 \) 不一定恰好是 \( \Delta t \) 的整数倍，代码采用 \( N_1 = \lceil T_1 / \Delta t \rceil \) 向上取整，再令 \( T_1 = N_1 \cdot \Delta t \) 进行精确调整：

\[
\Delta t = \frac{T}{N_t}, \quad N_1 = \left\lceil \frac{T_1}{\Delta t} \right\rceil, \quad T_1 \leftarrow N_1 \cdot \Delta t, \quad T_x = T + T_1
\]

代码中 `dt = T/Nt`，`N1 = int(math.ceil(T1/dt))`，`T1 = N1*dt`，`Tx = T + T1`，完全对应上式。

---

### 4.2 乘积矩形规则的离散方程

**对应代码（fouvol.py，第278-291行：solver==1，基本乘积矩形规则）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L278)

```python
# fouvol.py 第278-291行：乘积矩形规则（solver==1）——实现公式 (12)
# product rectangle rule
if solver == 1:
  # set up the denominator for product rectangle rule
  if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
  # begin time stepping
  for i in range(1,Nt+1):
    # get history, list comprehension seems around 5% faster
    hist = sum([varphi0/(alpha+1)*(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*uvals[j] for j in range(1,i)])
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)
```

**论文中的乘积矩形规则离散方程**

乘积矩形规则（Product Rectangle Rule）对 Volterra 积分方程 (10) 进行如下离散：寻找 \( U_n \approx u(t_n) \)，对 \( n = 1, 2, \ldots \) 逐步求解，满足：

\[
U_n + \sum_{j=1}^{n} U_j \int_{t_{j-1}}^{t_j} \phi(t_n - s)\, ds = f(t_n)
\]

其中被积函数 \( \phi(t_n - s) = \phi_0(t_n - s + \bar{t})^\alpha \) 在每段 \( [t_{j-1}, t_j] \) 上视 \( U_j \) 为常数（矩形近似），积分精确计算。代码中 `for i in range(1,Nt+1)` 逐步推进 \( n = 1, 2, \ldots, N_t \)，`hist` 积累历史项，`uvals[i]` 存储 \( U_n \)。

---

### 4.3 乘积积分的精确公式

**对应代码（fouvol.py，第287行：hist 的计算表达式）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L287)

```python
# fouvol.py 第287行：乘积矩形规则历史项的核心计算——对应积分精确计算公式
hist = sum([varphi0/(alpha+1)*(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*uvals[j] for j in range(1,i)])
```

**论文中的积分精确计算公式**

在乘积矩形规则中，利用 \( \alpha \neq -1 \) 的条件，对每段区间 \( [t_{j-1}, t_j] \) 上的幂律核进行精确积分：

\[
\phi_0 \int_{t_{j-1}}^{t_j} \left(\bar{t} + t_n - s\right)^{\alpha} ds = \frac{\phi_0}{\alpha+1} \left[ \left(\bar{t} + t_n - t_{j-1}\right)^{\alpha+1} - \left(\bar{t} + t_n - t_j\right)^{\alpha+1} \right]
\]

代入 \( t_n = n\Delta t \)，\( t_j = j\Delta t \)，\( t_{j-1} = (j-1)\Delta t \)：

\[
= \frac{\phi_0}{\alpha+1} \left[ \left(\bar{t} + (n-j+1)\Delta t\right)^{\alpha+1} - \left(\bar{t} + (n-j)\Delta t\right)^{\alpha+1} \right]
\]

代码直接实现此公式：
- `varphi0/(alpha+1)` 对应 \( \phi_0/(\alpha+1) \)
- `((i-j+1)*dt+tbar)**(alpha+1)` 对应 \( (\bar{t} + (n-j+1)\Delta t)^{\alpha+1} \)
- `((i-j)*dt+tbar)**(alpha+1)` 对应 \( (\bar{t} + (n-j)\Delta t)^{\alpha+1} \)
- 乘以 `uvals[j]`（即 \( U_j \)）后在 `j=1` 到 `i-1` 上求和

---

### 4.4 分母 d 的计算公式

**对应代码（fouvol.py，第279-283行：denom 的设置）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L279)

```python
# fouvol.py 第279-283行：分母 d 的计算，对应论文中 d 的定义
# set up the denominator for product rectangle rule
if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
```

**论文中分母 \( d \) 的定义公式**

在乘积矩形规则中，\( j = n \) 时的积分项移至方程左侧，形成分母：

\[
d := 1 + \int_{t_{n-1}}^{t_n} \phi(t_n - s)\, ds = 1 + \frac{\phi_0}{\alpha+1}\left[\left(\bar{t} + \Delta t\right)^{\alpha+1} - \bar{t}^{\alpha+1}\right]
\]

注意 \( d \) 与时间步 \( n \) 无关（因为 \( \phi \) 仅依赖 \( t_n - s \)，而此段长度固定为 \( \Delta t \)），所以可在时间循环外预计算。代码中：
- `tbar > 0` 分支：`denom = 1.0 + varphi0/(alpha+1)*((dt+tbar)**(alpha+1)-tbar**(alpha+1))`，对应 \( 1 + \frac{\phi_0}{\alpha+1}[(\bar{t}+\Delta t)^{\alpha+1} - \bar{t}^{\alpha+1}] \)
- `tbar == 0` 分支（弱奇异情形）：`denom = 1.0 + varphi0/(alpha+1)*(dt**(alpha+1))`，即 \( 1 + \frac{\phi_0}{\alpha+1}\Delta t^{\alpha+1} \)（因为 \( \bar{t} = 0 \) 时 \( \bar{t}^{\alpha+1} = 0 \)）

---

### 4.5 公式 (12)：基本乘积矩形规则的完整求解公式

**对应代码（fouvol.py，第288-290行：uvals[i] 的更新）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L288)

```python
# fouvol.py 第288-290行：乘积矩形规则的时间步进——直接实现公式 (12)
    hist = sum([varphi0/(alpha+1)*(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*uvals[j] for j in range(1,i)])
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)
```

**论文公式 (12)：基本乘积矩形规则**

将第4.3节的积分精确计算公式和第4.4节的分母 \( d \) 代入乘积矩形规则，得到完整的 \( U_n \) 求解公式：

\[
\begin{aligned}
U_n &= \frac{f(t_n) - \displaystyle\sum_{j=1}^{n-1} U_j \int_{t_{j-1}}^{t_j} \phi(t_n - s)\, ds}{1 + \displaystyle\int_{t_{n-1}}^{t_n} \phi(t_n - s)\, ds} \\[10pt]
&= \frac{f(t_n) - \displaystyle\sum_{j=1}^{n-1} \frac{\phi_0 U_j}{\alpha+1}\left[(\bar{t}+(n-j+1)\Delta t)^{\alpha+1} - (\bar{t}+(n-j)\Delta t)^{\alpha+1}\right]}{1 + \dfrac{\phi_0}{\alpha+1}\left[(\bar{t}+\Delta t)^{\alpha+1} - \bar{t}^{\alpha+1}\right]} \\[10pt]
&= d^{-1} f(t_n) - d^{-1} \sum_{j=1}^{n-1} \frac{\phi_0 U_j}{\alpha+1}\left[(\bar{t}+(n-j+1)\Delta t)^{\alpha+1} - (\bar{t}+(n-j)\Delta t)^{\alpha+1}\right] \quad (12)
\end{aligned}
\]

代码中 `uvals[i] = (ft(i*dt, varphi0, alpha) - hist) / denom` 精确实现此公式：
- `ft(i*dt, varphi0, alpha)` 计算 \( f(t_n) \)
- `hist` 是历史求和项
- `denom` 是分母 \( d \)
- 整体表达式即 \( d^{-1}f(t_n) - d^{-1}\cdot\text{hist} \)

---

### 4.6 乘积梯形规则（solver==2）

**对应代码（fouvol.py，第293-314行：solver==2，乘积梯形规则）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L293)

```python
# fouvol.py 第293-314行：乘积梯形规则（solver==2）——高阶基本规则
# product trapezoidal rule for all except first interval
elif solver == 2:
  # set up the denominator for product rectangle rule for first interval
  if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
  # product rectangle rule on first interval, convenience setting for u[0] (used in history sum)
  uvals[1] = ft(dt,varphi0,alpha) / denom
  uvals[0] = uvals[1]
  # set up the denominator for product trapezoidal rule
  if tbar > 0:
    denom = 1.0 + 0.5*varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + 0.5*varphi0/(alpha+1)*( dt**(alpha+1) )
  # begin time stepping
  for i in range(2,Nt+1):
    hist = sum([0.5*varphi0/(alpha+1)*(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*(uvals[j]+uvals[j-1]) \
                   for j in range(1,i)])
    hist = hist + 0.5*varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) ) * uvals[i-1]
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)
```

**乘积梯形规则的积分近似**

乘积梯形规则在每段 \( [t_{j-1}, t_j] \) 上用梯形近似 \( u(s) \approx \frac{1}{2}(U_j + U_{j-1}) \)，而对核 \( \phi \) 仍精确积分：

\[
\int_{t_{j-1}}^{t_j} \phi(t_n - s)\, u(s)\, ds \approx \frac{U_j + U_{j-1}}{2} \cdot \frac{\phi_0}{\alpha+1}\left[(\bar{t}+(n-j+1)\Delta t)^{\alpha+1} - (\bar{t}+(n-j)\Delta t)^{\alpha+1}\right]
\]

梯形规则的分母变为：

\[
d_{\text{trap}} := 1 + \frac{1}{2} \cdot \frac{\phi_0}{\alpha+1}\left[(\bar{t}+\Delta t)^{\alpha+1} - \bar{t}^{\alpha+1}\right]
\]

代码 `0.5*varphi0/(alpha+1)*...*(uvals[j]+uvals[j-1])` 对应上式（`0.5` 即梯形规则的 \( 1/2 \) 因子，`uvals[j]+uvals[j-1]` 对应 \( U_j + U_{j-1} \)）。第一步使用矩形规则初始化（`uvals[1] = ft(dt,...) / denom`），第二步起采用梯形规则。

---

### 4.7 公式 (13)-(14)：引入 Fourier 代理后的历史分割

**对应代码（fouvol.py，第317-353行：solver==3，矩形-Fourier 法完整实现）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L317)

```python
# fouvol.py 第317-353行：矩形-Fourier 法（solver==3），实现公式 (13)-(16)
# Fourier, with product rectangle rule
elif solver == 3:
  # set up the denominator for product rectangle rule
  if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
  # begin time stepping
  for i in range(1,Nt+1):
    # update Fourier history variables from previous time step
    if i > N1:
      for j in range(-Nlim,0):
        Rpart = phistarR[j+Nlim]; Ipart = phistarI[j+Nlim];
        costerm = math.cos(math.pi*j*dt/Tx); sinterm = math.sin(math.pi*j*dt/Tx)
        cosT1   = math.cos(math.pi*j*T1/Tx); sinT1   = math.sin(math.pi*j*T1/Tx)
        # complex attenuation
        phistarR[j+Nlim] = costerm*Rpart - sinterm*Ipart
        phistarI[j+Nlim] = sinterm*Rpart + costerm*Ipart
        # product integration improvement - the old way used the rectangle rule
        update = cosT1*sinterm-sinT1*(1.0-costerm)
        coeff = cnR[j+Nlim]*uvals[i-N1]*Tx/j/math.pi
        phistarR[j+Nlim] = phistarR[j+Nlim] + coeff*update
        update = cosT1*(1.0-costerm)+sinT1*sinterm
        phistarI[j+Nlim] = phistarI[j+Nlim] + coeff*update
      # the special case where j=0 - no 'attenuation' is necessary
      phistarR[0+Nlim] = phistarR[0+Nlim] + cnR[0+Nlim]*uvals[i-N1]*dt

    # quadrature over j={i-N1,...,i-2,i-1}
    jlo = int(max(1,i+1-N1))
    hist = varphi0/(alpha+1)*sum([(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*uvals[j] for j in range(jlo,i)])
    # Fourier contribution for far history
    if i > N1:
      hist = hist + phistarR[0+Nlim]
      hist = hist + 2*sum([phistarR[j+Nlim] for j in range(-Nlim,0)])
    # solve...
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)
```

**论文公式 (13)-(14)：历史分割与 Fourier 展开**

当 \( n > N_1 \) 时，将历史求和分成近端和远端两部分（公式13）：

\[
\begin{aligned}
& \left(1 + \int_{t_{n-1}}^{t_n} \phi(t_n - s)\, ds\right) U_n \\
= & \; f(t_n) - \sum_{j=n-N_1+1}^{n-1} U_j \int_{t_{j-1}}^{t_j} \phi(t_n - s)\, ds - \sum_{j=1}^{n-N_1} U_j \int_{t_{j-1}}^{t_j} \psi(t_n - s)\, ds \quad (13)
\end{aligned}
\]

将 \( \psi \) 用截断 Fourier 级数代理（公式14），其中 \( L \in \mathbb{N} \) 由 `Nlim` 指定：

\[
\approx f(t_n) - \sum_{j=n-N_1+1}^{n-1} U_j \int_{t_{j-1}}^{t_j} \phi(t_n - s)\, ds - \sum_{k=0}^{L} c_k \sum_{j=1}^{n-N_1} U_j \int_{t_{j-1}}^{t_j} e^{i\pi k (t_n - s)/T_x}\, ds \quad (14)
\]

代码中：
- 近端历史（`jlo` 到 `i-1` 的求和中 `hist`）对应公式(13)的第二项（\( \sum_{j=n-N_1+1}^{n-1} \)）
- 远端 Fourier 历史（`phistarR` 数组）对应公式(14)的最后一项

当 `i <= N1` 时，`if i > N1` 分支不执行，退化为基本矩形规则（公式12）。

---

### 4.8 \(\mathcal{H}_k(n)\) 的定义

**对应代码（fouvol.py，第273-275行：phistarR 和 phistarI 的初始化）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L273)

```python
# fouvol.py 第273-275行：复 Fourier 历史变量的初始化（H_k(N1) = 0）
# set up the complex Fourier history - whether or not it is needed
phistarR = np.zeros(2*Nlim+1) 
phistarI = np.zeros(2*Nlim+1) 
```

**\(\mathcal{H}_k(n)\) 的数学定义**

论文定义了辅助量 \( \mathcal{H}_k(n) \)，它是 Fourier 代理对远端历史积分的贡献：

\[
\mathcal{H}_k(n) := c_k \sum_{j=1}^{n-N_1} U_j \int_{t_{j-1}}^{t_j} e^{i\pi k (t_n - s)/T_x}\, ds
\]

以及

\[
\mathcal{H}_k(n-1) := c_k \sum_{j=1}^{(n-1)-N_1} U_j \int_{t_{j-1}}^{t_j} e^{i\pi k (t_{n-1} - s)/T_x}\, ds
\]

代码中 `phistarR[j+Nlim]` 存储 \( \mathcal{H}_k(n) \) 的**实部**，`phistarI[j+Nlim]` 存储**虚部**（索引 \( j \in [-L, L] \)，偏移量为 `Nlim`）。初始时 `phistarR = phistarI = 0`（即 \( \mathcal{H}_k(N_1) = 0 \)，因为在 \( n \leq N_1 \) 时尚无远端历史）。

---

### 4.9 公式 (15)：Fourier 历史变量的递归更新

**对应代码（fouvol.py，第325-341行：复指数衰减与乘积积分更新）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L325)

```python
# fouvol.py 第325-341行：递归更新 H_k(n)——实现论文公式 (15) 的核心
    if i > N1:
      for j in range(-Nlim,0):
        Rpart = phistarR[j+Nlim]; Ipart = phistarI[j+Nlim];
        costerm = math.cos(math.pi*j*dt/Tx); sinterm = math.sin(math.pi*j*dt/Tx)
        cosT1   = math.cos(math.pi*j*T1/Tx); sinT1   = math.sin(math.pi*j*T1/Tx)
        # complex attenuation
        phistarR[j+Nlim] = costerm*Rpart - sinterm*Ipart
        phistarI[j+Nlim] = sinterm*Rpart + costerm*Ipart
        # product integration improvement - the old way used the rectangle rule
        update = cosT1*sinterm-sinT1*(1.0-costerm)
        coeff = cnR[j+Nlim]*uvals[i-N1]*Tx/j/math.pi
        phistarR[j+Nlim] = phistarR[j+Nlim] + coeff*update
        update = cosT1*(1.0-costerm)+sinT1*sinterm
        phistarI[j+Nlim] = phistarI[j+Nlim] + coeff*update
      # the special case where j=0 - no 'attenuation' is necessary
      phistarR[0+Nlim] = phistarR[0+Nlim] + cnR[0+Nlim]*uvals[i-N1]*dt
```

**论文公式 (15)：递归更新公式**

公式 (15) 是本算法的核心递推公式。利用

\[
e^{i\pi k (t_n - s)/T_x} = e^{i\pi k \Delta t / T_x} \cdot e^{i\pi k (t_{n-1} - s)/T_x}
\]

可将 \( \mathcal{H}_k(n) \) 分解为新增项与旧历史衰减更新之和：

\[
\begin{aligned}
\mathcal{H}_k(n) &= c_k U_{n-N_1} \int_{t_{n-N_1-1}}^{t_{n-N_1}} e^{i\pi k (t_n - s)/T_x}\, ds \\
&\quad + e^{i\pi k \Delta t / T_x} \cdot c_k \sum_{j=1}^{(n-1)-N_1} U_j \int_{t_{j-1}}^{t_j} e^{i\pi k (t_{n-1} - s)/T_x}\, ds \\
&= c_k U_{n-N_1} \int_{t_{n-N_1-1}}^{t_{n-N_1}} e^{i\pi k (t_n - s)/T_x}\, ds + e^{i\pi k \Delta t / T_x} \cdot \mathcal{H}_k(n-1) \quad (15)
\end{aligned}
\]

**代码与公式的对应关系：**

**第一步——复指数衰减（attenuation）**，即乘以 \( e^{i\pi k \Delta t/T_x} \)：

\[
\mathcal{H}_k(n) \leftarrow e^{i\pi k \Delta t / T_x} \cdot \mathcal{H}_k(n-1)
\]

用复数乘法展开（令 \( \omega = e^{i\pi k\Delta t/T_x} = \cos(\pi k\Delta t/T_x) + i\sin(\pi k\Delta t/T_x) \)）：

\[
\begin{pmatrix}\text{Re}\,\mathcal{H}_k(n) \\ \text{Im}\,\mathcal{H}_k(n)\end{pmatrix} \leftarrow \begin{pmatrix}\cos\theta & -\sin\theta \\ \sin\theta & \cos\theta\end{pmatrix} \begin{pmatrix}\text{Re}\,\mathcal{H}_k(n-1) \\ \text{Im}\,\mathcal{H}_k(n-1)\end{pmatrix}
\]

其中 \( \theta = \pi k \Delta t / T_x \)。代码中：

```python
costerm = math.cos(math.pi*j*dt/Tx)   # cos(pi k dt / Tx)
sinterm = math.sin(math.pi*j*dt/Tx)   # sin(pi k dt / Tx)
# 实部更新：Re = cos*Re_old - sin*Im_old
phistarR[j+Nlim] = costerm*Rpart - sinterm*Ipart
# 虚部更新：Im = sin*Re_old + cos*Im_old
phistarI[j+Nlim] = sinterm*Rpart + costerm*Ipart
```

**第二步——新增项的乘积积分**，即 \( c_k U_{n-N_1} \int_{t_{n-N_1-1}}^{t_{n-N_1}} e^{i\pi k(t_n - s)/T_x} ds \)：

对区间 \( [t_{n-N_1-1}, t_{n-N_1}] \) 精确积分（令 \( s_0 = t_{n-N_1-1} \)，\( s_1 = t_{n-N_1} \)，则 \( t_n - s_0 = (N_1+1)\Delta t \)，\( t_n - s_1 = N_1\Delta t = T_1 \)）：

\[
\int_{t_{n-N_1-1}}^{t_{n-N_1}} e^{i\pi k(t_n-s)/T_x} ds = \frac{T_x}{i\pi k}\left[e^{i\pi k(N_1+1)\Delta t/T_x} - e^{i\pi k T_1/T_x}\right]
\]

展开（用 \( \phi_1 = \pi k T_1 / T_x \) 和 \( \theta = \pi k \Delta t / T_x \)）：

\[
= \frac{T_x}{\pi k}\left[(\cos\phi_1\cdot\sin\theta - \sin\phi_1\cdot(1-\cos\theta)) + i(\cos\phi_1\cdot(1-\cos\theta) + \sin\phi_1\cdot\sin\theta)\right]
\]

代码中：

```python
cosT1   = math.cos(math.pi*j*T1/Tx)   # cos(pi k T1 / Tx)
sinT1   = math.sin(math.pi*j*T1/Tx)   # sin(pi k T1 / Tx)
# 实部更新分量
update = cosT1*sinterm - sinT1*(1.0-costerm)
coeff = cnR[j+Nlim]*uvals[i-N1]*Tx/j/math.pi   # c_k * U_{n-N1} * Tx / (k*pi)
phistarR[j+Nlim] = phistarR[j+Nlim] + coeff*update
# 虚部更新分量
update = cosT1*(1.0-costerm) + sinT1*sinterm
phistarI[j+Nlim] = phistarI[j+Nlim] + coeff*update
```

**特殊情况 \( k = 0 \)** 时，\( e^{i\pi \cdot 0 \cdot (t_n-s)/T_x} = 1 \)，积分退化为区间长度 \( \Delta t \)，无需衰减：

\[
\mathcal{H}_0(n) = \mathcal{H}_0(n-1) + c_0 \cdot U_{n-N_1} \cdot \Delta t
\]

代码中 `phistarR[0+Nlim] = phistarR[0+Nlim] + cnR[0+Nlim]*uvals[i-N1]*dt` 直接实现此公式。

---

### 4.10 公式 (16)：带 Fourier 代理的最终求解公式

**对应代码（fouvol.py，第343-351行：hist 与 Fourier 贡献之和）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L343)

```python
# fouvol.py 第343-351行：组合近端历史与 Fourier 远端历史，实现公式 (16)
    # quadrature over j={i-N1,...,i-2,i-1}
    jlo = int(max(1,i+1-N1))
    hist = varphi0/(alpha+1)*sum([(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*uvals[j] for j in range(jlo,i)])
    # Fourier contribution for far history
    if i > N1:
      hist = hist + phistarR[0+Nlim]
      hist = hist + 2*sum([phistarR[j+Nlim] for j in range(-Nlim,0)])
    # solve...
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)
```

**论文公式 (16)：带 Fourier 代理的完整求解公式**

综合公式 (14) 和 (15)，最终得到带 Fourier 代理的时间步进公式：

\[
U_n = d^{-1} f(t_n) - d^{-1} \sum_{j=n-N_1+1}^{n-1} U_j \int_{t_{j-1}}^{t_j} \phi(t_n - s)\, ds - d^{-1} \sum_{k=0}^{L} \mathcal{H}_k(n) \quad (16)
\]

代码的对应关系：
- `hist` 的第一部分（`jlo` 到 `i-1` 的求和）对应 \( \sum_{j=n-N_1+1}^{n-1} U_j \int_{t_{j-1}}^{t_j} \phi(t_n-s)ds \)（公式16第二项分子）
- `phistarR[0+Nlim]` 对应 \( \mathcal{H}_0(n) \)（\( k=0 \) 项，已乘以 \( c_0 \)）
- `2*sum([phistarR[j+Nlim] for j in range(-Nlim,0)])` 对应 \( \sum_{k=1}^{L}\mathcal{H}_k(n) + \sum_{k=-L}^{-1}\mathcal{H}_k(n) = 2\sum_{k=1}^{L}\text{Re}\,\mathcal{H}_k(n) \)（利用 \( \mathcal{H}_{-k} = \overline{\mathcal{H}_k} \) 的共轭对称性，实部之和等于两倍）
- `uvals[i] = (ft(...) - hist) / denom` 实现 \( U_n = d^{-1}(f(t_n) - \text{近端历史} - \text{远端Fourier历史}) \)

---

### 4.11 梯形-Fourier 规则（solver==4）

**对应代码（fouvol.py，第356-402行：solver==4，梯形-Fourier 法）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L356)

```python
# fouvol.py 第356-402行：梯形-Fourier 规则（solver==4）
# Fourier, with product trapezoidal rule
elif solver == 4:
  # set up the denominator for product rectangle rule for first interval
  if tbar > 0:
    denom = 1.0 + varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + varphi0/(alpha+1)*( dt**(alpha+1) )
  # product rectangle rule on first interval, convenience setting for u[0] (used in history sum)
  uvals[1] = ft(dt,varphi0,alpha) / denom
  uvals[0] = uvals[1]
  # set up the denominator for product trapezoidal rule
  if tbar > 0:
    denom = 1.0 + 0.5*varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )
  else:
    denom = 1.0 + 0.5*varphi0/(alpha+1)*( dt**(alpha+1) )
  # begin time stepping
  for i in range(2,Nt+1):
    # update Fourier history variables from previous time step
    if i > N1:
      for j in range(-Nlim,0):
        Rpart = phistarR[j+Nlim]; Ipart = phistarI[j+Nlim];
        costerm = math.cos(math.pi*j*dt/Tx); sinterm = math.sin(math.pi*j*dt/Tx)
        cosT1   = math.cos(math.pi*j*T1/Tx); sinT1   = math.sin(math.pi*j*T1/Tx)
        # complex attenuation
        phistarR[j+Nlim] = costerm*Rpart - sinterm*Ipart
        phistarI[j+Nlim] = sinterm*Rpart + costerm*Ipart
        # product integration
        update = cosT1*sinterm-sinT1*(1.0-costerm)
        coeff = cnR[j+Nlim]*(uvals[i-N1]+uvals[i-1-N1])/2.0*Tx/j/math.pi
        phistarR[j+Nlim] = phistarR[j+Nlim] + coeff*update
        update = cosT1*(1.0-costerm)+sinT1*sinterm
        phistarI[j+Nlim] = phistarI[j+Nlim] + coeff*update
      # the special case where j=0 - no 'attenuation' is necessary
      phistarR[0+Nlim] = phistarR[0+Nlim] + cnR[0+Nlim]*(uvals[i-N1]+uvals[i-1-N1])/2.0*dt

    # quadrature over j={i-N1,...,i-2,i-1}
    jlo = int(max(1,i+1-N1))
    hist = sum([0.5*varphi0/(alpha+1)*(((i-j+1)*dt+tbar)**(alpha+1)-((i-j)*dt+tbar)**(alpha+1))*(uvals[j]+uvals[j-1]) \
           for j in range(jlo,i)])
    hist = hist + 0.5*varphi0/(alpha+1)*( (dt+tbar)**(alpha+1)-tbar**(alpha+1) )*uvals[i-1]
    # Fourier contribution for far history
    if i > N1:
      hist = hist + phistarR[0+Nlim]
      hist = hist + 2*sum([phistarR[j+Nlim] for j in range(-Nlim,0)])
    # solve...
    uvals[i] = ( ft(i*dt,varphi0,alpha) - hist ) / denom
    if cheat and i<Nt/4:
      uvals[i] = ut(i*dt,alpha)
```

**梯形-Fourier 规则的递归更新（公式15的梯形版本）**

梯形-Fourier 规则将近端历史改用梯形规则，而远端 Fourier 代理部分也采用梯形近似 \( U_{n-N_1} \to \frac{1}{2}(U_{n-N_1} + U_{n-1-N_1}) \)：

对 \( k \neq 0 \)，公式 (15) 的梯形版本为：

\[
\mathcal{H}_k(n) = c_k \cdot \frac{U_{n-N_1} + U_{n-1-N_1}}{2} \int_{t_{n-N_1-1}}^{t_{n-N_1}} e^{i\pi k(t_n-s)/T_x} ds + e^{i\pi k\Delta t/T_x} \cdot \mathcal{H}_k(n-1)
\]

对 \( k = 0 \)：

\[
\mathcal{H}_0(n) = c_0 \cdot \frac{U_{n-N_1} + U_{n-1-N_1}}{2} \cdot \Delta t + \mathcal{H}_0(n-1)
\]

代码中 `coeff = cnR[j+Nlim]*(uvals[i-N1]+uvals[i-1-N1])/2.0*Tx/j/math.pi` 将矩形规则的 `uvals[i-N1]` 替换为梯形平均 `(uvals[i-N1]+uvals[i-1-N1])/2.0`，其余结构与 solver==3 完全一致。

---

### 4.12 时间步进算法总结

**对应代码（fouvol.py，第270-275行：初始化）**

[点击查看源代码位置](https://github.com/xingzhiyuan1229/-Fourier-series-recursion-for-problems-/blob/copilot/add-numerical-method-formula-interpretation/fouvol.py#L270)

```python
# fouvol.py 第270-275行：解向量与 Fourier 历史变量的初始化
# begin solution process, store solution here - avoid solving for u(0)
uvals = np.zeros(Nt+1)  # don't need the ,1) here - it makes abs() return [xxx] not xxx
cheat = 0
# set up the complex Fourier history - whether or not it is needed
phistarR = np.zeros(2*Nlim+1) 
phistarI = np.zeros(2*Nlim+1) 
```

**论文中的时间步进算法（结合公式12、15、16的完整步骤）**

论文第4节最后给出完整的时间步进算法：

1. 使用公式 (12) 计算 \( U_n \)（对 \( n = 1, 2, \ldots, N_1 \)）：

\[
U_n = d^{-1} f(t_n) - d^{-1} \sum_{j=1}^{n-1} \frac{\phi_0 U_j}{\alpha+1}\left[(\bar{t}+(n-j+1)\Delta t)^{\alpha+1} - (\bar{t}+(n-j)\Delta t)^{\alpha+1}\right]
\]

2. 初始化每个 \( \mathcal{H}_k(N_1) = 0 \)（代码：`phistarR = phistarI = np.zeros(2*Nlim+1)`）

3. 对 \( n = N_1+1, N_1+2, \ldots, N_t \)：
   - 使用公式 (15) 更新每个 \( \mathcal{H}_k(n) \)（代码：`phistarR`, `phistarI` 的衰减和新增更新）
   - 使用公式 (16) 计算 \( U_n \)（代码：`uvals[i] = (ft(...) - hist) / denom`）

代码中的 `if i > N1` 分支自然实现了步骤1到步骤3的切换，从基本矩形规则无缝过渡到 Fourier 代理算法。`uvals = np.zeros(Nt+1)` 初始化解向量（含 \( U_0 = 0 \)），`phistarR` 和 `phistarI` 分别存储 \( \mathcal{H}_k(n) \) 的实部和虚部，索引方式为 `[k + Nlim]`（支持 \( k \in [-L, L] \)）。

---

## 完整总结

本文档（第一部分与第二部分）覆盖了论文翻译.md 中第2节至第4节出现的全部公式，并给出了每个公式在代码库中最底层的对应计算实现，具体对应关系汇总如下：

| 论文公式 | 对应代码文件 | 对应函数/代码位置 | 核心实现 |
|---------|------------|-----------------|---------|
| 公式(10)：Volterra积分方程 | `fouvol.py` | 第87-99行（参数定义）| `varphi0`, `alpha`, `tbar` 等参数 |
| 精确解 \( u(t) = t^{-\alpha} \) | `fouvol.py` | 第29-30行 `ut()` | `return t**(-alpha)` |
| 公式(11)：右端函数 \( f(t) \) | `fouvol.py` | 第33-34行 `ft()` | `t**(-alpha) + varphi0*pi*(-alpha)*t/sin(-pi*alpha)` |
| 积分分割（第3节） | `fouvol.py` | 第343-351行（solver==3）| `hist` 的近端+远端分割 |
| Hermite插值边界条件 | `fouker.py` | 第282-331行 `get_FourierCoefficientsNew()` | `da`, `db` 导数计算 |
| 分差矩阵构造 | `fouker.py` | 第141-166行 `get_DDHermite()` | `DDmat[k,k+n] = da[n]/nfac` |
| Hermite插值求值 | `fouker.py` | 第169-175行 `get_DDHermiteValue()` | Newton插值公式 |
| Fourier系数 \( c_k \)（复指数形式） | `fouker.py` | 第333-352行 `get_FourierCoefficientsNew()` | `cnR[n+Nlim] = 1.0/Tx*realpart[0]` |
| Fourier系数（余弦形式） | `fouker.py` | 第336-350行 | `quad(... cos(pi*n*t/Tx) ...)` |
| 分母 \( d \)（公式12的分母） | `fouvol.py` | 第279-283行（solver==1）| `denom = 1.0 + varphi0/(alpha+1)*(...)` |
| 公式(12)：基本矩形规则 | `fouvol.py` | 第285-290行（solver==1）| `uvals[i] = (ft(...) - hist) / denom` |
| 公式(13)-(14)：Fourier分割 | `fouvol.py` | 第317-351行（solver==3）| `if i > N1:` 分支整体 |
| \( \mathcal{H}_k(n) \) 初始化 | `fouvol.py` | 第273-275行 | `phistarR = phistarI = np.zeros(...)` |
| 公式(15)：递归更新 | `fouvol.py` | 第325-341行（solver==3）| 复指数衰减 + 乘积积分新增项 |
| 公式(16)：最终求解 | `fouvol.py` | 第343-351行（solver==3）| `hist + 2*sum(phistarR...)` |
| 梯形-Fourier（公式15变体） | `fouvol.py` | 第356-402行（solver==4）| `(uvals[i-N1]+uvals[i-1-N1])/2.0` |
