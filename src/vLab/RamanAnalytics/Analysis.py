import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn import model_selection
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import Ridge, RidgeCV, Lasso, LassoCV, ElasticNet, ElasticNetCV, LinearRegression


class RR:
    """
    Ridge Regression:         
    Minimizes the objective function:

    .. math:: ||y - Xw||^2_2 + \\alpha * ||w||^2_2

    This model solves a regression model where the loss function is
    the linear least squares function with regularization is given by
    the l2-norm. 
    """

    def __init__(self):
        pass

    def analysis(self, X_train, X_test, y_train, y_test, Metabolite, def_alpha=None):
        """
        :param X_train: input data at training set;
        :param X_test: input data at test set;
        :param y_train: output at training set;
        :param y_test:  output at test set;
        :param Metabolite: the name of output;
        :param def_alpha: Regularization strength; must be a positive float.
        """
        if def_alpha == None:
            # Cross-validation
            ridgecv = RidgeCV(scoring='neg_mean_squared_error')
            ridgecv.fit(X_train, y_train)
            def_alpha = ridgecv.alpha_

        ridge = Ridge(alpha=def_alpha)
        ridge.fit(X_train, y_train)
        mean_squared_error(y_test, ridge.predict(X_test))

        ridge_coef = ridge.coef_

        plt.figure()
        g = plt.scatter(pd.to_numeric(X_train.columns), ridge_coef)
        g.axes.set_title('Regression Coefficients (RR)')
        g.axes.set_xlabel('Raman Shift (cm-1)')
        g.axes.set_ylabel('Coefficient')
        plt.show()

        predicted_RR_train = ridge.predict(X_train)

        plt.figure()
        g = plt.scatter(y_train, predicted_RR_train)
        g.axes.set_title(str(Metabolite) + ' Concentration Train (RR)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        # plt.plot(predicted_value_train,predicted_value_train, 'r-' )
        plt.show()

        predicted_RR_test = ridge.predict(X_test)

        plt.figure()
        g = plt.scatter(y_test, predicted_RR_test)
        g.axes.set_title(str(Metabolite) + ' Concentration Test (RR)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        # plt.plot([0,1],[0,np.max(g.axes.get_xlim(),g.axes.get_ylim()),predicted_value_test)], 'r-' )
        plt.show()


class LA:
    """
    Lasso: 
    Minimizes the objective function:
        
    .. math:: ||y - Xw||^2_2 + \\alpha * ||w||_1
    
    This model solves a regression model where the loss function is
    the linear least squares function with L1 prior as regularizer (aka the Lasso).    
    """

    def __init__(self):
        pass

    def analysis(self, X_train, X_test, y_train, y_test, Metabolite, def_alpha=None):
        """
        :param X_train: input data at training set;
        :param X_test: input data at test set;
        :param y_train: output at training set;
        :param y_test:  output at test set;
        :param Metabolite: the name of output;
        :param def_alpha: Regularization strength; must be a positive float.
        """
        if def_alpha is None:
            # Cross-validation
            lassoCV = LassoCV(max_iter=10000)
            lassoCV.fit(X_train, y_train)
            def_alpha = lassoCV.alpha_

        lasso = Lasso(alpha=def_alpha, normalize=True)
        lasso.fit(X_train, y_train)
        mean_squared_error(y_test, lasso.predict(X_test))

        lasso_coef = lasso.coef_

        plt.figure()
        g = plt.scatter(pd.to_numeric(X_train.columns), lasso_coef)
        g.axes.set_title('Regression Coefficients (LA)')
        g.axes.set_xlabel('Raman Shift (cm-1)')
        g.axes.set_ylabel('Coefficient')
        plt.show()

        predicted_LA_train = lasso.predict(X_train)

        plt.figure()
        g = plt.scatter(y_train, predicted_LA_train)
        g.axes.set_title(str(Metabolite) + ' Concentration Train (LA)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        plt.show()

        predicted_LA_test = lasso.predict(X_test)

        plt.figure()
        g = plt.scatter(y_test, predicted_LA_test)
        g.axes.set_title(str(Metabolite) + ' Concentration Test (LA)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        plt.show()


class EN:
    """
    Elastic Net:

    Minimizes the objective function:

    .. math:: ||y - Xw||^2_2 + \\alpha * r * ||w||_1 + 0.5 * \\alpha * (1 - r) * ||w||^2_2

    The parameter r is the L1 ratio. r = 1 is the lasso penalty and r = 0 is the ridge regression penalty.
    """

    def __init__(self):
        pass

    def analysis(self, X_train, X_test, y_train, y_test, Metabolite, def_alpha=None):
        """
        :param X_train: input data at training set;
        :param X_test: input data at test set;
        :param y_train: output at training set;
        :param y_test:  output at test set;
        :param Metabolite: the name of output;
        :param def_alpha: Regularization strength; must be a positive float.
        """
        if def_alpha == None:
            # Cross-validation
            elasticCV = ElasticNetCV(max_iter=10000)
            elasticCV.fit(X_train, y_train)
            def_alpha = elasticCV.alpha_

        elastic = ElasticNet(alpha=def_alpha, normalize=True)
        elastic.fit(X_train, y_train)
        mean_squared_error(y_test, elastic.predict(X_test))

        elastic_coef = elastic.coef_

        plt.figure()
        g = plt.scatter(pd.to_numeric(X_train.columns), elastic_coef)
        g.axes.set_title('Regression Coefficients (EN)')
        g.axes.set_xlabel('Raman Shift (cm-1)')
        g.axes.set_ylabel('Coefficient')
        plt.show()

        predicted_EN_train = elastic.predict(X_train)

        plt.figure()
        g = plt.scatter(y_train, predicted_EN_train)
        g.axes.set_title(str(Metabolite) + ' Concentration Train (EN)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        plt.show()

        predicted_EN_test = elastic.predict(X_test)

        plt.figure()
        g = plt.scatter(y_test, predicted_EN_test)
        g.axes.set_title(str(Metabolite) + ' Concentration Test (EN)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        plt.show()


class PCR:
    """
    Principal component analysis (PCA) with Linear Regression.
    
    PCA:  
    Linear dimensionality reduction using Singular Value Decomposition of the
    data to project it to a lower dimensional space. The input data is centered
    but not scaled for each feature before applying the SVD.
    
    Linear Regression:   
    LinearRegression fits a linear model with coefficients 
    to minimize the residual sum of squares between the observed targets in
    the dataset, and the targets predicted by the linear approximation.
    
    """

    def __init__(self):
        pass

    def analysis(self, X_train, X_test, y_train, y_test, Metabolite, def_comp=None):
        """
        :param X_train: input data at training set;
        :param X_test: input data at test set;
        :param y_train: output at training set;
        :param y_test:  output at test set;
        :param Metabolite: the name of output;
        :param def_comp: number of components to keep
        """
        pca = PCA()

        # Scale the data
        X_reduced_train = pca.fit_transform(scale(X_train))
        n = len(X_reduced_train)

        if def_comp == None:
            # 10-fold CV, with shuffle
            kf_10 = model_selection.KFold(n_splits=10, shuffle=True, random_state=1)

            mse = []
            regr = LinearRegression()
            # Calculate MSE with only the intercept (no principal components in regression)
            score = -1 * model_selection.cross_val_score(regr, np.ones((n, 1)), y_train, cv=kf_10,
                                                         scoring='neg_mean_squared_error').mean()
            mse.append(score)

            # Calculate MSE using CV for the 20 principle components, adding one component at the time.
            for i in np.arange(1, 21):
                score = -1 * model_selection.cross_val_score(regr, X_reduced_train[:, :i], y_train, cv=kf_10,
                                                             scoring='neg_mean_squared_error').mean()
                mse.append(score)

            # plt.plot(np.array(mse), '-v')
            # plt.xlabel('Number of principal components in regression')
            # plt.ylabel('MSE')
            # plt.title(str(Metabolite))
            # plt.xlim(xmin=-1);

            def_comp = mse.index(np.min(mse)) + 1

        # Train regression model on training data 
        regr = LinearRegression()
        regr.fit(X_reduced_train[:, :(def_comp)], y_train)

        pcr_coef = np.matmul(pca.components_.T[:, 0:def_comp], regr.coef_.T)

        plt.figure()
        g = plt.scatter(pd.to_numeric(X_train.columns), pcr_coef)
        g.axes.set_title('Regression Coefficients (PCR)')
        g.axes.set_xlabel('Raman Shift (cm-1)')
        g.axes.set_ylabel('Coefficient')
        plt.show()

        # test = np.matmul(scale(X_train), np.matmul(pca.components_.T[:,0:PCA_component], regr.coef_.T)) + np.mean(y_train)[0]
        # np.matmul(scale(X_train), pca.components_.T[:,0:PCA_component])

        predicted_PCR_train = regr.predict(X_reduced_train[:, :(def_comp)])

        plt.figure()
        g = plt.scatter(y_train, predicted_PCR_train)
        g.axes.set_title(str(Metabolite) + ' Concentration Train (PCR)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        plt.show()

        # Prediction with test data
        X_reduced_test = pca.transform(scale(X_test))[:, :(def_comp)]
        predicted_PCR_test = regr.predict(X_reduced_test)
        mean_squared_error(y_test, predicted_PCR_test)

        plt.figure()
        g = plt.scatter(y_test, predicted_PCR_test)
        g.axes.set_title(str(Metabolite) + ' Concentration Test (PCR)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        plt.show()


class PLS:
    """
    Partial least squares regression (PLS):
    find the multidimensional direction in the input space that explains 
    the maximum multidimensional variance direction in the output space.
    """

    def __init__(self):
        pass

    def analysis(self, X_train, X_test, y_train, y_test, Metabolite, def_comp=None):
        """
        :param X_train: input data at training set;
        :param X_test: input data at test set;
        :param y_train: output at training set;
        :param y_test:  output at test set;
        :param Metabolite: the name of output;
        :param def_comp: number of components to keep
        """
        # n = len(X_train)

        if def_comp is None:
            # 10-fold CV, with shuffle
            kf_10 = model_selection.KFold(n_splits=10, shuffle=True, random_state=1)

            mse = []

            for i in np.arange(1, 20):
                pls = PLSRegression(n_components=i)
                score = -1 * model_selection.cross_val_score(pls, scale(X_train), y_train, cv=kf_10,
                                                             scoring='neg_mean_squared_error').mean()
                mse.append(score)

            # # Plot results
            # plt.plot(np.arange(1, 20), np.array(mse), '-v')
            # plt.xlabel('Number of principal components in regression')
            # plt.ylabel('MSE')
            # plt.title(str(Metabolite))
            # plt.xlim(xmin=-1)

            def_comp = mse.index(np.min(mse)) + 1

        pls = PLSRegression(n_components=def_comp)
        pls.fit(scale(X_train), y_train)

        pls_coef = pls.coef_

        plt.figure()
        g = plt.scatter(pd.to_numeric(X_train.columns), pls_coef)
        g.axes.set_title('Regression Coefficients (PLS)')
        g.axes.set_xlabel('Raman Shift (cm-1)')
        g.axes.set_ylabel('Coefficient')
        plt.show()

        predicted_PLS_train = pls.predict(scale(X_train))

        plt.figure()
        g = plt.scatter(y_train, predicted_PLS_train)
        g.axes.set_title(str(Metabolite) + ' Concentration Train (PLS)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        plt.show()

        # Prediction with test data
        predicted_PLS_test = pls.predict(scale(X_test))
        mean_squared_error(y_test, pls.predict(scale(X_test)))

        plt.figure()
        g = plt.scatter(y_test, predicted_PLS_test)
        g.axes.set_title(str(Metabolite) + ' Concentration Test (PLS)')
        g.axes.set_xlabel('True Values (g/L)')
        g.axes.set_ylabel('Predictions (g/L)')
        g.axes.axis('equal')
        g.axes.axis('square')
        g.axes.axline([0, 0], [1, 1], color='r')
        plt.show()
