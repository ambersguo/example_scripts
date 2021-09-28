import sys
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model.logistic import LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import roc_curve, auc, confusion_matrix
import matplotlib.pyplot as plt

with open('Logistic_regression_model_evaluations.txt','w') as file2out:
    sys.stdout = file2out

    df = pd.read_csv('./inclone_nonclone_diseases_for_logistic_regression_for_sklearn.txt', delimiter='\t')
    X_train_raw, X_test_raw, y_train, y_test = train_test_split(df['diseases'], df['inclone'], random_state=11)
    vectorizer = TfidfVectorizer()
    X_train = vectorizer.fit_transform(X_train_raw)
    X_test = vectorizer.transform(X_test_raw)
    #solver: for small datasets, 'liblinear' is a good choice, whereas 'sag' and'saga' are faster for large ones
    classifier = LogisticRegression(solver='liblinear')
    classifier.fit(X_train, y_train)
    #LogisticRegression.score() predicts and scores labels for a test set using accuracy
    scores = cross_val_score(classifier, X_train, y_train, cv=5)
    print('Accuracy: %s' % np.mean(scores))

    precisions = cross_val_score(classifier, X_train, y_train, cv=5, scoring='precision')
    print('Precision: %s' % np.mean(precisions))
    recalls = cross_val_score(classifier, X_train, y_train, cv=5, scoring='recall')
    print('Recall: %s' % np.mean(recalls))

    f1s = cross_val_score(classifier, X_train, y_train, cv=5, scoring='f1')
    print('F1 score: %s' % np.mean(f1s))
    
    #create a performance matrix to evaluate logistic classifier
    y_pred = classifier.predict(X_test)
    confusion_matrix = confusion_matrix(y_test, y_pred)
    print('confusion matrix:')
    print(confusion_matrix)
    plt.matshow(confusion_matrix)
    plt.title('Confusion matrix')
    plt.colorbar()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.savefig('binary_classification_with_logistic_regression_sklearn_confusion_matrix.pdf')
    plt.close()

    #AUC is the area under the ROC curve; it reduces the ROC curve to a single value,
    #which represents the expected performance of the classifier. The dashed line in the
    #following figure is for a classifier that predicts classes randomly; it has an AUC of
    #0.5. The solid curve is for a classifier that outperforms random guessing
    predictions = classifier.predict_proba(X_test)
    false_positive_rate, recall, thresholds = roc_curve(y_test, predictions[:, 1])
    roc_auc = auc(false_positive_rate, recall)
    plt.title('Receiver Operating Characteristic')
    plt.plot(false_positive_rate, recall, 'b', label='AUC = %0.2f' % roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.ylabel('Recall')
    plt.xlabel('Fall-out')
    plt.savefig('binary_classification_with_logistic_regression_sklearn_ROC.pdf')
    plt.close()

    def show_most_informative_features(vectorizer, clf, n=20):
        feature_names = vectorizer.get_feature_names()
        coefs_with_fns = sorted(zip(clf.coef_[0], feature_names))
        top = zip(coefs_with_fns[:n], coefs_with_fns[:-(n + 1):-1])
        print('Coefficient of feature not in clone\tTop feature not in clone\tCoefficient of feature in clone\tTop feature in clone')
        for (coef_1, fn_1), (coef_2, fn_2) in top:
            print ('%.4f\t%-15s\t%.4f\t%-15s' % (coef_1, fn_1, coef_2, fn_2))

    show_most_informative_features(vectorizer, classifier)
