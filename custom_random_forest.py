import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier
from multiprocessing import Pool


def fit_tree(X, y, i, max_features, max_depth, random_state):
    """
    Fits a DecisionTreeClassifier on a bootstrap sample of the data.

    Parameters:
    - X : array-like or sparse matrix, shape (n_samples, n_features)
        The training input samples.
    - y : array-like, shape (n_samples,)
        The target values.
    - i : int
        Index of the tree being fitted.
    - max_features : int
        The number of features to consider when looking for the best split.
    - max_depth : int
        The maximum depth of the tree.
    - random_state : int
        Random seed for reproducibility.

    Returns:
    - tree : DecisionTreeClassifier
        The fitted decision tree.
    - feat_ids : array-like
        Indices of selected features for the tree.
    """
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
    """
    Predicts class probabilities for the input samples using a given decision tree.

    Parameters:
    - X : array-like or sparse matrix, shape (n_samples, n_features)
        The input samples.
    - feat_ids : array-like
        Indices of selected features for the tree.
    - tree : DecisionTreeClassifier
        The trained decision tree.

    Returns:
    - proba : array-like, shape (n_samples, n_classes)
        Class probabilities of the input samples.
    """
    proba = tree.predict_proba(X[:, feat_ids])
    return proba


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
            self, n_estimators=10, max_depth=None, max_features=None, random_state=42
    ):
        """
        Initializes the RandomForestClassifierCustom instance.

        Parameters:
        - n_estimators : int, default=10
            The number of trees in the forest.
        - max_depth : int, default=None
            The maximum depth of the trees. If None, nodes are expanded until
            all leaves are pure or until all leaves contain less than
            min_samples_split samples.
        - max_features : int, default=None
            The number of features to consider when looking for the best split.
        - random_state : int, default=42
            Random seed for reproducibility.
        """
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []


    def fit(self, X, y, n_jobs=1):
        """
        Fits the RandomForestClassifierCustom to the training data.

        Parameters:
        - X : array-like or sparse matrix, shape (n_samples, n_features)
            The training input samples.
        - y : array-like, shape (n_samples,)
            The target values.
        - n_jobs : int, default=1
            The number of jobs to run in parallel for fitting the trees.

        Returns:
        - self : RandomForestClassifierCustom
            Returns the instance itself.
        """
        self.classes_ = sorted(np.unique(y))
        
        with Pool(n_jobs) as pool:
            poolup = list(pool.starmap(fit_tree, [(X, y, i, self.max_features, self.max_depth, self.random_state) for i in range(self.n_estimators)]))
            for i in poolup:
                self.trees.append(i[0])
                self.feat_ids_by_tree.append(i[1])
            
        return self

    def predict_proba(self, X, n_jobs=1):
        """
        Predicts class probabilities for the input samples.

        Parameters:
        - X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples.
        - n_jobs : int, default=1
            The number of jobs to run in parallel.

        Returns:
        - avg_proba : array-like, shape (n_samples, n_classes)
            Average class probabilities across all trees.
        """    
        with Pool(n_jobs) as pool:
            probas = pool.starmap(predict_proba_tree, [(X, self.feat_ids_by_tree[i], self.trees[i]) for i in range(len(self.trees))])

        avg_proba = np.mean(probas, axis=0)
        return avg_proba

    def predict(self, X, n_jobs=1):
        """
        Predicts class labels for the input samples.

        Parameters:
        - X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples.
        - n_jobs : int, default=1
            The number of jobs to run in parallel.

        Returns:
        - predictions : array-like, shape (n_samples,)
            Predicted class labels.
        """
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
