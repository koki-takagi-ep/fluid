# 差分公式のテイラー展開による導出

本ドキュメントでは、有限差分法で使用する差分公式をテイラー展開から導出する。

## 目次

1. [テイラー展開の基本](#1-テイラー展開の基本)
2. [1階微分の差分公式](#2-1階微分の差分公式)
3. [2階微分の差分公式](#3-2階微分の差分公式)
4. [風上差分の導出](#4-風上差分の導出)

---

## 1. テイラー展開の基本

関数 $f(x)$ の点 $x_0$ まわりのテイラー展開：

$$
f(x_0 + h) = f(x_0) + h f'(x_0) + \frac{h^2}{2!} f''(x_0) + \frac{h^3}{3!} f'''(x_0) + \frac{h^4}{4!} f''''(x_0) + O(h^5)
$$

$$
f(x_0 - h) = f(x_0) - h f'(x_0) + \frac{h^2}{2!} f''(x_0) - \frac{h^3}{3!} f'''(x_0) + \frac{h^4}{4!} f''''(x_0) + O(h^5)
$$

以下、$f_i = f(x_i)$、$h = \Delta x$ と記す。

---

## 2. 1階微分の差分公式

### 2.1 前進差分（Forward Difference）

$f(x_0 + h)$ のテイラー展開を $f'(x_0)$ について解く：

$$
f(x_0 + h) = f(x_0) + h f'(x_0) + \frac{h^2}{2} f''(x_0) + O(h^3)
$$

$$
f'(x_0) = \frac{f(x_0 + h) - f(x_0)}{h} - \frac{h}{2} f''(x_0) + O(h^2)
$$

したがって：

$$
\boxed{\frac{\partial f}{\partial x}\bigg|_i \approx \frac{f_{i+1} - f_i}{\Delta x} + O(\Delta x)}
$$

**精度: 1次**

### 2.2 後退差分（Backward Difference）

$f(x_0 - h)$ のテイラー展開を $f'(x_0)$ について解く：

$$
f(x_0 - h) = f(x_0) - h f'(x_0) + \frac{h^2}{2} f''(x_0) + O(h^3)
$$

$$
f'(x_0) = \frac{f(x_0) - f(x_0 - h)}{h} + \frac{h}{2} f''(x_0) + O(h^2)
$$

したがって：

$$
\boxed{\frac{\partial f}{\partial x}\bigg|_i \approx \frac{f_i - f_{i-1}}{\Delta x} + O(\Delta x)}
$$

**精度: 1次**

### 2.3 中心差分（Central Difference）

前進展開と後退展開を**引き算**する：

$$
f(x_0 + h) - f(x_0 - h) = 2h f'(x_0) + \frac{2h^3}{3!} f'''(x_0) + O(h^5)
$$

$$
f'(x_0) = \frac{f(x_0 + h) - f(x_0 - h)}{2h} - \frac{h^2}{6} f'''(x_0) + O(h^4)
$$

したがって：

$$
\boxed{\frac{\partial f}{\partial x}\bigg|_i \approx \frac{f_{i+1} - f_{i-1}}{2\Delta x} + O(\Delta x^2)}
$$

**精度: 2次**（1次の誤差項が相殺）

---

## 3. 2階微分の差分公式

### 3.1 2次中心差分

前進展開と後退展開を**足し算**する：

$$
f(x_0 + h) + f(x_0 - h) = 2f(x_0) + h^2 f''(x_0) + \frac{h^4}{12} f''''(x_0) + O(h^6)
$$

$f''(x_0)$ について解く：

$$
f''(x_0) = \frac{f(x_0 + h) - 2f(x_0) + f(x_0 - h)}{h^2} - \frac{h^2}{12} f''''(x_0) + O(h^4)
$$

したがって：

$$
\boxed{\frac{\partial^2 f}{\partial x^2}\bigg|_i \approx \frac{f_{i+1} - 2f_i + f_{i-1}}{\Delta x^2} + O(\Delta x^2)}
$$

**精度: 2次**

### 3.2 ラプラシアン（2次元）

2次元の場合、各方向に上記の公式を適用：

$$
\nabla^2 f = \frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2}
$$

$$
\boxed{\nabla^2 f \bigg|_{i,j} \approx \frac{f_{i+1,j} - 2f_{i,j} + f_{i-1,j}}{\Delta x^2} + \frac{f_{i,j+1} - 2f_{i,j} + f_{i,j-1}}{\Delta y^2} + O(\Delta x^2, \Delta y^2)}
$$

**精度: 2次**（各方向で2次精度）

---

## 4. 風上差分の導出

### 4.1 移流方程式と特性線

移流方程式 $\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = 0$ において、情報は速度 $c$ の符号に応じた方向から伝播する。

- $c > 0$: 情報は左（$x$ 負方向）から来る → 後退差分を使用
- $c < 0$: 情報は右（$x$ 正方向）から来る → 前進差分を使用

### 4.2 1次風上差分

$$
c \frac{\partial u}{\partial x} \approx c \times \begin{cases}
\frac{u_i - u_{i-1}}{\Delta x} & (c \geq 0) \\[8pt]
\frac{u_{i+1} - u_i}{\Delta x} & (c < 0)
\end{cases}
$$

**精度: 1次**

### 4.3 数値拡散

1次風上差分のテイラー展開：

$$
\frac{u_i - u_{i-1}}{\Delta x} = \frac{\partial u}{\partial x} - \frac{\Delta x}{2} \frac{\partial^2 u}{\partial x^2} + O(\Delta x^2)
$$

移流方程式に代入すると：

$$
\frac{\partial u}{\partial t} + c \frac{\partial u}{\partial x} = \frac{c \Delta x}{2} \frac{\partial^2 u}{\partial x^2} + O(\Delta x^2)
$$

右辺の $\frac{c \Delta x}{2} \frac{\partial^2 u}{\partial x^2}$ は**数値拡散**（numerical diffusion）であり、物理的な拡散ではない人工的な拡散を発生させる。

---

## 差分公式のまとめ

| 差分公式 | 式 | 精度 |
|---------|-----|------|
| 前進差分 | $\frac{f_{i+1} - f_i}{\Delta x}$ | $O(\Delta x)$ |
| 後退差分 | $\frac{f_i - f_{i-1}}{\Delta x}$ | $O(\Delta x)$ |
| 中心差分（1階） | $\frac{f_{i+1} - f_{i-1}}{2\Delta x}$ | $O(\Delta x^2)$ |
| 中心差分（2階） | $\frac{f_{i+1} - 2f_i + f_{i-1}}{\Delta x^2}$ | $O(\Delta x^2)$ |

---

## 参考文献

- LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and Partial Differential Equations*. SIAM.
- Ferziger, J. H., & Perić, M. (2002). *Computational Methods for Fluid Dynamics* (3rd ed.). Springer.
