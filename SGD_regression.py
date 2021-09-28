import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from sklearn.linear_model import SGDRegressor
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

with open('SGD_model_evaluations.txt','w') as file2out:
    sys.stdout = file2out

    #using pandas to load txt data to a dataframe
    df = pd.read_csv('clone_size_for_SGD.txt', delimiter = '\t')
    #calculates summary statistics for each column of the dataframe
    print(df.describe())
    #calculates a pairwise correlation matrix to show correlation among all features
    #USAGE::DataFrame.corr(self, method='pearson', min_periods=1), method : {‘pearson’, ‘kendall’, ‘spearman’} or callable
    print(df.corr())

    features = df.drop('clone_size', axis = 1)
    X_train, X_test, y_train, y_test = train_test_split(features, df.clone_size)
    #scale the features using StandardScaler
    X_scaler = StandardScaler()
    y_scaler = StandardScaler()
    y_train = np.array(y_train).reshape(-1,1)
    y_test = np.array(y_test).reshape(-1,1)
    X_train = X_scaler.fit_transform(X_train)
    y_train = y_scaler.fit_transform(y_train)
    X_test = X_scaler.transform(X_test)
    y_test = y_scaler.transform(y_test)
    #stochastic gradient descent(SGD) to optimize the residual sum of squares(cost function)
    regressor = SGDRegressor(loss='squared_loss',max_iter=1000,tol=0.001)
    scores = cross_val_score(regressor, X_train, y_train.ravel(), cv=5)
    print('Cross validation r-squared scores: ' + str(scores) +'\n')
    print('Average cross validation r-squared score: ' + str(np.mean(scores)) +'\n')
    regressor.fit(X_train, y_train.ravel())
    print('Test set r-squared score: ' + str(regressor.score(X_test, y_test.ravel())) +'\n')
    print('Best Intercept and Coef: ' + str(regressor.intercept_) + str(regressor.coef_) +'\n')
