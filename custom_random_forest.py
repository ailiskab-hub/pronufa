import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier
from multiprocessing import Pool


def fit_tree(X, y, i, max_features, max_depth, random_state):
    np.random.seed(random_state + i)
    n_samples, n_features = X.shape

    feat_ids = np.random.choice(n_features, max_features, replace=False)

    bs_ind = np.random.choice(n_samples, n_samples, replace=True)
    X_bs = X[bs_ind, :]
    y_bs = y[bs_ind]

    tree = DecisionTreeClassifier(max_depth=max_depth, random_state=random_state)

    tree.fit(X_bs[:, feat_ids], y_bs)
    return tree, feat_ids


def predict_proba_tree(X, feat_ids, tree):
    proba = tree.predict_proba(X[:, feat_ids])
    return proba


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
            self, n_estimators=10, max_depth=None, max_features=None, random_state=42
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []


    def fit(self, X, y, n_jobs=1):
        self.classes_ = sorted(np.unique(y))
        
        with Pool(n_jobs) as pool:
            poolup = list(pool.starmap(fit_tree, [(X, y, i, self.max_features, self.max_depth, self.random_state) for i in range(self.n_estimators)]))
            for i in poolup:
                self.trees.append(i[0])
                self.feat_ids_by_tree.append(i[1])
            
        return self

    def predict_proba(self, X, n_jobs=1):
    
        with Pool(n_jobs) as pool:
            probas = pool.starmap(predict_proba_tree, [(X, self.feat_ids_by_tree[i], self.trees[i]) for i in range(len(self.trees))])

        avg_proba = np.mean(probas, axis=0)
        return avg_proba

    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
