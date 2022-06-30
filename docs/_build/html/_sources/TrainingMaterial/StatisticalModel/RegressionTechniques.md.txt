This section introduces several regression methods, including Least Absolute Shrinkage and Selection Operator (lasso), Ridge Regression (RR), Elastic Net (EN),  Partial Least Squares (PLS), and Principal Component Regression (PCR). After finish learning this section, the students should have ability to choose appropriate methods, tuning the hyperparameters, train the model and interpretate results.

# Regression Techniques
The regression models considered are linear models with the general form:

$$
y = \beta_0 x_0 + \beta_1 x_1 + \cdots + \beta_{n-1} x_{n-1} + \epsilon,
$$

where the $y$,  $x_i$ and $\beta_i$ represent output, $i$-th inputs  and $i$-th regression coefficients; $x_0 = 1$; $n$ is the number of inputs and regression coefficients; and $\epsilon$ denotes the error term. Without loss of generality, this model can be rewritten in matrix form as

$$
\bf y = \bf X \bf \beta + \bf{\epsilon},
$$

where the input data matrix $\bf{X} \in \mathcal{R}^{m \times n}$, the vector of regression coefficients $\bf{\beta} \in \mathcal{R}^{n \times 1}$, $\bf{y} \in \mathcal{R}^{m \times 1}$, and $m$ is the number of observations. The elements of $\bf{X}$ can be any mix of raw input data and transformations of the raw data (aka
features). The errors $\bf \epsilon$ are assumed to be homoscedastic, have zero mean, and are uncorrelated.


Model building involves determining the vector $\bf \beta$ from the data $\bf X$ and $\bf y$ that minimizes the error $\bf \epsilon$ concerning a defined measure of the error.

## Ordinary Least Squares Regression
The key idea of OLS regression is to find a solution $\bf \beta$  that minimizes the L2-norm of the model error $\bf \epsilon$ (Strang, 2016),

$$
min_\bf \beta || \bf y - \bf X \bf \beta||_2^2
$$

Practically, this optimization can be solved and leads to the analytical solution

$$
\bf{\beta} = \bf X^T (\bf X \bf X ^T)^{-1} \bf y
$$

## Lasso
OLS estimates may have low bias but have very large variance for many real-world data analytics problems, resulting in low prediction accuracy on unseen data. Lasso is a strategy for addressing this problem by adding the L1-norm of the weights as a penalty to the least-squares objective (Tibshirani, 1996),

$$
min_\bf \beta || \bf y - \bf X \bf \beta||_2^2 + \lambda ||\bf{\beta}||_1,
$$

For positive values for $\lambda$, lasso builds a model in which some weights are zero, i.e., $\beta_i = 0$, due to the structure of the optimization. The larger the value of $\lambda$, the more weights are forced to zero. Therefore, Lasso can be used to selects a subset of variables to be used in the regression and can potentially lead to increased interpretability of results by removing measurements that are not needed in the model prediction.


## Ridge Regression
The underlying motivation for Ridge Regression (RR) is identical to lasso. The difference between the methods is that Ridge Regression (RR) adds a L2-norm penalty to the objective function,

$$
min_\bf \beta || \bf y - \bf X \bf \beta||_2^2 + \lambda ||\bf{\beta}||_2^2,
$$

Setting the derivative of the objective function to zero leads to the closed-form solution [Hoerl and Kennard, 1970]

$$
\bf{\beta} =  (\bf X \bf X ^T + \lambda \bf{I})^{-1} \bf X^T \bf y
$$

As in lasso, the L2-norm penalty on$\lambda$ can be rewritten as a constraint. Due to the shape of the constraint in the resulting optimization, the weights $\beta_i$ in RR will never reach zero although they can be arbitrarily small [James et al., 2021]. Similarly to lasso, RR improves the predictive accuracy of the model by introducing a bias that reduces variance in the estimated parameters [Zou and Hastie, 2005]. However, models produced by RR can be challenging to interpret, since all model inputs are retained in the model even for high values of $\lambda$.

## Elastic Net
The Elastic Net (EN) is a combination of lasso and RR, 

$$
min_\bf \beta || \bf y - \bf X \bf \beta||_2^2 + \lambda P(\bf{\beta}),
$$

$$
P(\bf{\beta}) = \frac{1-\alpha}{2} ||\bf{\beta}||_2^2 + \alpha ||\bf{\beta}||_1
$$

for $\alpha \in (0, 1)$ and non-negative values of $\lambda$. The L1-norm provides the ability to set unimportant parameters to zero whereas the L2-norm improves robustness in the selection of which parameters to retain in the model. Lasso and RR are limiting cases of the EN, for $\alpha$ equal to 1 and 0, respectively.

## PCA and PCR
Principal Component Analysis (PCA) is an approach to manipulating the data matrix $\bf{X} \in \mathcal{R}^{m \times n}$. The idea of PCA is to find a lower dimensional representation of $\bf X$ while conserving the variations in the data. From an optimization perspective, PCA solves

$$
max_{||w||_2=1} ||\bf X w||^2_2,
$$

for each principal component, with subsequent principal components also required to be orthogonal to the previous principal components. PCA finds a weighted combination of the mean subtracted columns of $\bf X$ retaining maximal variance. For the first principal component, the constrained optimization is equivalent to the unconstrained optimization: 

$$
w_1 = argmax_w \frac{w^T X^T X w}{w^Tw},
$$ 

The optimization is solved by $w$ being the eigenvector corresponding to the largest eigenvalue of the positive semidefinite matrix $\bf{X}^T \bf{X}$.

The complete projection step of PCA can be written as 

$$
\bf{T} = \bf{X W}
$$ 

where $\bf{T} \in \mathcal{R}^{m \times n}$ is the PCA score matrix, $\bf{X} \in \mathcal{R}^{m \times n}$ is the data matrix, and $\bf{W} \in \mathcal{R}^{n \times n}$ is a coefficient matrix. The columns of $\bf{W}$ are called loadings and correspond to the eigenvectors of $\bf{X^T X}$ in descending order.

The above equation describes a linear transformation. The loadings form an orthogonal basis, which results in the columns of $\bf T$ being decorrelated (Bengio et al., 2013). The first $l < n$ components of $\bf T$ form the $l$-dimensional representation of $\bf{X}$ preserving most variance and are denoted by $\bf{T}_l$.

By combining the OLS with PCA leads to Principal Component Regression (PCR):

$$
\bf{\beta}_l = \bf T_l^T (\bf T_l \bf T_l ^T)^{-1} \bf y,
$$

where $l$ is the number of principal components. Consequently, the regression coefficients in the original space are constructed
by: 

$$
\bf{\beta} = \bf{W}_l \bf{\beta}_l,
$$

where $\bf{W}_l$ denotes the first $l$ columns of the matrix $\bf{W}$.

## PLS
Partial Least Squares (PLS) aims to find lower dimensional representations of $\bf X$ and $\bf Y$ and is not restricted to scalar objectives. While PCA performs dimensionality reduction in an unsupervised way, PLS incorporates information about the target $\bf Y$ in the dimensionality reduction scheme. In general, the governing equations for PLS are

$$
\bf X = \bf{TP}^T + \bf{E},
$$

$$
\bf Y = \bf{UQ}^T + \bf{F},
$$

where the matrix $\bf{X} \in \mathcal{R}^{m \times n}$ is the data matrix, and $\bf{Y} \in \mathcal{R}^{m \times p}$ is the matrix of responses, $\bf{T} \in \mathcal{R}^{m \times l}$ is the score matrix, $\bf{P} \in \mathcal{R}^{n \times l}$ is the loading matrix, and $\bf{E} \in \mathcal{R}^{m \times n}$ is the residual matrix corresponding to $\bf{X}$. Similarly for $\bf{Y}$, the matrix $\bf{U} \in \mathcal{R}^{m \times l}$ is its score matrix, $\bf{Q} \in \mathcal{R}^{p \times l}$ is its loading matrix, and $\bf{F} \in \mathcal{R}^{m \times p}$ is its residual matrix.

Similarly to PCA, PLS performs a linear transformation on the input $\bf{X}, with 

$$
\bf{T} = \bf{XW}
$$

From an optimization perspective, PLS maximizes the sample covariance between the X scores and the responses (Wold et al., 2001, Boulesteix and Strimmer, 2006). The calculation of the first component can be written in the form of the unconstrained optimization

$$
\bf{w}_1 = argmax_w \frac{\bf{w}^T \bf{X}^T y y \bf{X} w}{w^Tw},
$$ 

Then the same regression can be applied as PCR.

## Reference

1. Gilbert Strang. Introduction to Linear Algebra. Cambridge Press, Wellesley, Massachusetts, fifth edition, 2016. ISBN
9780980232776.
2. Robert Tibshirani. Regression shrinkage and selection via the lasso. Journal of the Royal Statistical Society: Series B (Methodological),
58(1):267–288, 1996. ISSN 0035-9246. doi:10.1111/j.2517-6161.1996.tb02080.x.
3. Hui Zou and Trevor Hastie. Regularization and variable selection via the elastic net. Journal of the Royal Statistical Society. Series
B: Statistical Methodology, 67(2):301–320, 2005. ISSN 13697412. doi:10.1111/j.1467-9868.2005.00503.x.
4. Arthur E. Hoerl and Robert W. Kennard. Ridge regression: Biased estimation for nonorthogonal problems. Technometrics, 12(1):
55–67, 1970. ISSN 15372723. doi:10.1080/00401706.1970.10488634.
5. Gareth James, Daniela Witten, Trevor Hastie, and Robert Tibshirani. Linear Model Selection and Regularization. Springer, New
York, 2021. ISBN 978-1-0716-1417-4. doi:10.1007/978-1-0716-1418-1.
6. Yoshua Bengio, Aaron Courville, and Pascal Vincent. Representation learning: A review and new perspectives. IEEE Transactions
on Pattern Analysis and Machine Intelligence, 35(8):1798–1828, 2013. doi:10.1109/TPAMI.2013.50.
7. Schaeffer, J., & Braatz, R. (2022). Latent Variable Method Demonstrator--Software for Understanding Multivariate Data Analytics Algorithms. arXiv preprint arXiv:2205.08132.
