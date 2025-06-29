# 最速降线的推导——变分法与 Euler-Lagrange 方程
> 简便起见, 本文在解决最速降线问题时均仅考虑了导数为 $0$, 没有证明极值性乃至最值性. 

## 最速降线问题
1696 年, 约翰·伯努利在写给他哥哥雅克布·伯努利的一封公开信中提出一个问题:
> 一个均匀重力场中有 $A, B$ 两点, $A$ 点处有一质点, 仅在重力作用下沿一曲线轨道滑向 $B$ 点, 问什么样的曲线能使得滑落时间最短.

他也解答了这道题, 这便是 `最速降线问题`. \
为了理解 _Bernoulli_ 给出的证明, 我们回顾几个物理学的定律.
> _*Law 1.* (Fermat 原理)_&emsp;光线传播的路径是需时最少的路径. \
> _*Law 2.* (折射定律)_&emsp;光在由一种介质传入另一介质后, 入射角 $i$ 和折射角 $\gamma$ 的正弦之比为一定值, 即折射率.

有这两个引理之后, 我们引入一引理.
> _*Lemma.*_&emsp;平面上有 $A, B$ 两点分布一条直线 $l$ 两侧, 给定常数 $v_1, v_2$，设点 $C\in l$, 过 $C$ 作 $l'\perp l$, 则 $\dfrac{|AC|}{v_1}+\dfrac{|CB|}{v_2}$ 取得极值当 $AC$ 与 $l'$ 夹角 $\phi_1$ 的正弦与 $BC$ 与 $l'$ 夹角 $\phi_2$ 的正弦之比为 $\dfrac{v_1}{v_2}$.

_*Proof.*_&emsp;不妨在平面直角坐标系上表示, 等比缩放坐标轴并平移、旋转使 $A(0, a), B(1, -b), l:y=0\  (a, b\in \mathbb{R}_+)$, 设 $C(0, x)\ (x\in(0, 1))$.

![fig1](https://cdn.luogu.com.cn/upload/image_hosting/df9u6myr.png)

设 $s(x)=v_1^{-1}|AC|+v_2^{-1}|CB|=v_1^{-1}\sqrt{x^2+a^2}+v_2^{-1}\sqrt{(1-x)^2+b^2}$, 则

$$
\begin{aligned}
s'(x)&=-v_1^{-1}\cdot\dfrac{1}{2}\cdot\dfrac{2x}{\sqrt{x^2+a^2}}-v_2^{-1}\cdot\dfrac{1}{2}\cdot\dfrac{2x-2}{\sqrt{(1-x)^2+b^2}} \\
&=v_2^{-1}\cdot\dfrac{1-x}{\sqrt{(1-x)^2+b^2}}-v_1^{-1}\cdot\dfrac{x}{\sqrt{x^2+a^2}} \\
&=\dfrac{\sin{\phi_2}}{v_2}-\dfrac{\sin{\phi_1}}{v_1}.
\end{aligned}
$$
取得极值时 $s'(x)=0$, 引理得证.&emsp;$\blacksquare$

回到最速降线问题, 如图所示. 当质点下滑时, 有动能定理
$$
mgy=\dfrac12mv^2 \Rightarrow v=\sqrt{2gy}.
$$

![fig2](https://cdn.luogu.com.cn/upload/image_hosting/ppsjbe4z.png)

为了保证滑落时间最短, 应该使得每一段微小竖直位移 $dy$ 中均有 $\dfrac{v_y}{\sin\theta_y}$ 为一定值, 其中 $\theta_y$ 为切线与 $y$ 轴的夹角, 即
$$
\dfrac{v_y}{\sin\theta_y}=C_{onstant}, \\
\Rightarrow\dfrac{v_y}{(1+y'^2)^{-\frac12}}=C_{onstant}, \\
\Rightarrow\sqrt{gy(1+y'^2)}=C_{onstant}, \\
\Rightarrow y(1+y'^2)=C_{onstant}.
$$
经变量分离解该微分方程的特解, 为一圆摆线.&emsp;$\blacksquare$

## 变分法与 _Euler-Lagrange_ 方程
_Bernoulli_ 的方法从微小位移的角度解决了这一问题. 如果从宏观角度描述这一问题呢?\
应用上文的动能定理, 重新表述最速降线问题如下:
> (_Bernoulli_ 最速降线问题的数学表述) 平面上两点 $A(0,0), B(1, b)$, 考虑函数系 $F=\{f|f(0)=0, f(1)=b\}$, 定义 
> $$
> \forall f\in F, T(f, f')=\int_0^1 \dfrac{\sqrt{1+f'(x)^2}}{\sqrt{f(x)}}\mathrm{d}x.
> $$
> 求使 $T(f, f')$ 取到极值的 $f$.


为了继续解决这一问题, 我们引入 _Euler-Lagrange_ 方程:
> _*Theorem.* (Euler-Lagrange)方程_&emsp;对于泛函 $\mathcal{L}(x, f(x), f'(x))$ 满足 $f(x_1), f(x_2)$ 均为定值, $\mathcal{F}(f)=\int_{x_1}^{x_2}\mathcal{L}(x, f, f')\mathrm{d}x$ 取得极值时一定有 $\dfrac{\partial\mathcal{L}}{\partial f}=\dfrac{\mathrm{d}}{\mathrm{d}x}(\dfrac{\partial\mathcal{L}}{\partial f'})$.

我们仍然可以仿照上述的费马原理证明法, 那便是尝试让 $T$ 对 $f$ 求导并取 $0$. 然而我们并不能直接对函数 "求导", 我们需要对 $f$ 进行一些变换, 找到一种可以求导的式来表示 $f$. 不妨回顾导数: 导数的本质就是研究对一个初始值和一个微小变化量, 函数值的变化情况.\
不妨考虑任意一个 $f$ 作为初始值, 设为 $\hat{f}$, 随后对于任意一个扰动函数 $\eta$ 满足 $\eta(x_1)=\eta(x_2)=0$, 将扰动函数以一系数 $\epsilon$ 叠加在 $f$ 上. 这样我们可以用对 $\epsilon$ 求导代替对 $f$ 的 "求导".\
_思考_: 如果初始函数 $\hat{f}$ 恰恰就是所求的最速降线, 那么求导有什么特点呢?\
对于任意的 $\hat{f}$, 求导为 $0$ 时的 $\hat{f}+\epsilon\eta$ 就是对应函数系得极值点, 那么最速降线就是整个 $F$ 的极值点, 它一定有 $\epsilon=0$ 是所有扰动函数 $\eta$ 对应的 $\hat{f}+\epsilon\eta$ 的极值点.\
那么
$$
\begin{aligned}
\dfrac{\mathrm{d}\mathcal{F}}{\mathrm{d}\epsilon}&=\int_{x_1}^{x_2}\dfrac{\mathrm{d}\mathcal{L}(x, \hat{f}+\epsilon\eta, \hat{f}'+\epsilon\eta')}{\mathrm{d}\epsilon}\mathrm{d}x \\
&=\int_{x_1}^{x_2}(\dfrac{\partial\mathcal{L}}{\partial f}\cdot\dfrac{\mathrm{d}(\hat{f}+\epsilon\eta)}{\mathrm{d}\epsilon}+\dfrac{\partial\mathcal{L}}{\partial f'}\cdot\dfrac{\mathrm{d}(\hat{f}'+\epsilon\eta')}{\mathrm{d}\epsilon})\mathrm{d}x \\
&=\int_{x_1}^{x_2}(\dfrac{\partial\mathcal{L}}{\partial f}\cdot\eta+\dfrac{\partial\mathcal{L}}{\partial f'}\cdot\eta')\mathrm{d}x \\
&=\int_{x_1}^{x_2}\dfrac{\partial\mathcal{L}}{\partial f}\cdot\eta\mathrm{d}x+\int_{x_1}^{x_2}\dfrac{\partial\mathcal{L}}{\partial f'}\cdot\eta'\mathrm{d}x \\
&=\int_{x_1}^{x_2}\dfrac{\partial\mathcal{L}}{\partial f}\cdot\eta\mathrm{d}x+(\dfrac{\partial\mathcal{L}}{\partial f'}\cdot\eta\bigg|_{x_1}^{x_2})-\int_{x_1}^{x_2}(\dfrac{\partial\mathcal{L}}{\partial f'})'\cdot\eta\mathrm{d}x,
\end{aligned}
$$
由于 $\eta(x_1)=\eta(x_2)=0$ 我们有 $\dfrac{\partial\mathcal{L}}{\partial f'}\cdot\eta\bigg|_{x_1}^{x_2}=0$, 进而
$$
\begin{aligned}
\dfrac{\mathrm{d}\mathcal{F}}{\mathrm{d}\epsilon}&=\int_{x_1}^{x_2}\dfrac{\partial\mathcal{L}}{\partial f}\cdot\eta\mathrm{d}x-\int_{x_1}^{x_2}(\dfrac{\partial\mathcal{L}}{\partial f'})'\cdot\eta\mathrm{d}x \\
&=\int_{x_1}^{x_2}\dfrac{\partial\mathcal{L}}{\partial f}\cdot\eta-(\dfrac{\partial\mathcal{L}}{\partial f'})'\cdot\eta\mathrm{d}x \\
&=\int_{x_1}^{x_2}(\dfrac{\partial\mathcal{L}}{\partial f}-(\dfrac{\partial\mathcal{L}}{\partial f'})')\cdot\eta\mathrm{d}x.
\end{aligned}
$$
我们要求 $\forall \eta, \dfrac{\mathrm{d}\mathcal{F}}{\mathrm{d}\epsilon}=0$, 于是有
$$
\dfrac{\partial\mathcal{L}}{\partial f}-(\dfrac{\partial\mathcal{L}}{\partial f'})'=0.
$$
此即为 _Euler-Lagrange_ 方程.&emsp;$\blacksquare$

## 应用 _Euler-Lagrange_ 方程解决最速降线问题
结论当然可以直接带入 _Euler-Lagrange_ 方程求解, 本文侧重求解 $\mathcal{L}$ 不显式含 $x$ 时的通用结论.
> _*Lemma.*_&emsp;$\mathcal{L}$ 不显式含 $x$ 时, _Euler-Lagrange_ 方程即 $\dfrac{\partial \mathcal{L}}{\partial f'}\cdot f'-\mathcal{L}=C_{onstant}$.

_*Proof.*_&emsp;两侧求导, 即证 $\dfrac{\mathrm{d}}{\mathrm{d}x}(\dfrac{\partial \mathcal{L}}{\partial f'}\cdot f'-\mathcal{L})=0$.\
即 $\dfrac{\mathrm{d}}{\mathrm{d}x}(\dfrac{\partial \mathcal{L}}{\partial f'}\cdot f')=\dfrac{\mathrm{d}\mathcal{L}}{\mathrm{d}x}$.
由 _Euler-Lagrange_  方程我们有
$$
\begin{aligned}
LHS&=f'\cdot\dfrac{\mathrm{d}}{\mathrm{d}x}(\dfrac{\partial \mathcal{L}}{\partial f'})+f''\cdot\dfrac{\partial\mathcal{L}}{\partial f'} \\
&=f'\cdot\dfrac{\partial \mathcal{L}}{\partial f}+f''\cdot\dfrac{\partial\mathcal{L}}{\partial f'}.
\end{aligned}
$$
另一方面,
$$
\begin{aligned}
RHS &= \dfrac{\partial\mathcal{L}}{\partial x}+\dfrac{\partial\mathcal{L}}{\partial f}\cdot f+\dfrac{\partial\mathcal{L}}{\partial f'}\cdot f' \\
&= f\cdot\dfrac{\partial\mathcal{L}}{\partial f}+f'\cdot\dfrac{\partial\mathcal{L}}{\partial f'} = LHS.\quad\blacksquare
\end{aligned}
$$

回归最速降线问题, 带入 _Euler-Lagrange_ 方程, 即
$$
f'\cdot\dfrac{\partial\frac{\sqrt{1+f'^2}}{\sqrt{2gf}}}{\partial f'}-\frac{\sqrt{1+f'^2}}{\sqrt{2gf}}=C_{onstant},
$$
化简得
$$
(1+f'^2)f=C_{onstant},
$$
与上文的微分方程相同, 其解即为圆摆线.&emsp;$\blacksquare$

## 附: 圆摆线微分方程的求解
$$
\begin{aligned}
(1+f'^2)f=C &\Rightarrow f'=\sqrt{\dfrac Cf-1} \\
&\Rightarrow \sqrt{\dfrac{f}{C-f}}\mathrm{d}f=\mathrm{d}x \\
&\Rightarrow \int\sqrt{\dfrac{f}{C-f}}\mathrm{d}f=\int\mathrm{d}x.
\end{aligned}
$$
令 $f=\frac C2 (1-\cos\theta)$ 有 $\mathrm{d}f=\frac C2\sin\theta\mathrm{d}\theta$.
$$
\begin{aligned}
\int\sqrt{\dfrac{f}{C-f}}\mathrm{d}f &= \frac C2 \int\sqrt{\dfrac{1-\cos\theta}{1+\cos\theta}}\sin\theta\mathrm{d}\theta \\
&= \frac C2 \int\sqrt{\dfrac{(1-\cos\theta)(1-\cos^2\theta)}{1+\cos\theta}}\mathrm{d}\theta \\
&=\frac C2 \int(1-\cos\theta)\mathrm{d}\theta \\
&= \frac C2(\theta+\sin\theta)=\int\mathrm{d}x=x+C_2.
\end{aligned}
$$
带入条件 $A(0, 0)$ 有
$$
\begin{cases}
x=\dfrac C2(\theta+\sin\theta) \\ \\
y=\dfrac C2(1-\cos\theta)
\end{cases}
$$
即圆摆线方程.&emsp;$\blacksquare$
