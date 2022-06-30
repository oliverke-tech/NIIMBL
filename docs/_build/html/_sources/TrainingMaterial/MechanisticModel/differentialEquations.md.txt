## Types of Differential Equations

Differential equations can be divided into several types. Apart from describing the properties of the equation itself, these classes of differential equations can help inform the choice of approach to a solution. Commonly used distinctions include whether the equation is

- ordinary or partial
- linear or nonlinear
- homogeneous or heterogeneous

There are many other properties and subclasses of differential equations which can be very useful in specific contexts, for examble, stable or non-stable solution of differential equations.

Classically, the first-order ordinary differential equation can be written as

$$
\begin{align}
\frac{dx}{dt}=f(t,x)
\end{align}
$$

### Equation Order

Differential equations are described by their order, determined by the term with the highest derivvative. An equation containing only first derivatives is a *first-order differential equation*, an equation containing the second derivative is a *second-order differential equation*, and so on.

### Homogeneous or Heterogeneous

- **Homogeneous differential equations** involve only $x$  and its derivatives of $t$, i.e.
  $\frac{dx}{dt}=f(x)$.
- **Nonhomogeneous differential equations** are the same as homogeneous differential equations, except they can have terms involving only *x* on the right side, as in this equation:$\frac{dx}{dt}=f(t,x)$.

### Ordinary or Partial DE

- *Ordinary differential equations* or (ODE) are equations where the derivatives are taken with respect to only one variable. That is, there is only one independent variable.  For example, simple cell culture equations can be written as

  $$
  \begin{align*}
  \frac{dx}{dt}&=\mu x\\
  \frac{dp}{dt}&=q x
  \end{align*}
  $$

  where $X$ and $P$ are viable cell density and product concentration, $\mu$ is cell growth rate, $q$ is specific production rate.

- *Partial differential equations* or (PDE) are equations that depend on partial derivatives of several variables. That is, there are several independent variables.
    - Transport equation: $\frac{\partial y}{dt}+c\frac{\partial y}{\partial x}=0$
    - Heat equation: $\frac{\partial y}{dt}=\frac{\partial^2 y}{\partial x^2}$


### Linear or Nonlinear

- **linear differential equation** is a differential equation that is defined by a linear polynomial
  in the unknown function and its derivatives, that is an equation of the form

  $$
  a_0x+a_1x^\prime+a_2 x^{\prime\prime}+\ldots + a_n x^{(n)}=b(t)
  $$

  where $x^{(n)}$ represent the $n$-th order derivative of $x$ with respect to $t$.

- A **non-linear differential equation** is a differential equation that is not a linear equation in the unknown function and its derivatives.