from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, RepeatedKFold, learning_curve, validation_curve, KFold, \
    StratifiedKFold, RepeatedStratifiedKFold, TimeSeriesSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import log_loss, make_scorer
from sklearn.metrics import accuracy_score, confusion_matrix, recall_score, precision_score, f1_score
import joblib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_learning_curve(train_sizes, train_scores, test_scores, nTrees, nDepth, nLeaf, alpha=0.25):

    # Training set mean and standard deviation
    train_mean = np.mean(train_scores, axis=1)
    train_std = np.std(train_scores, axis=1)

    # Test set mean and standard deviation
    test_mean = np.mean(test_scores, axis=1)
    test_std = np.std(test_scores, axis=1)

    # Plot
    plt.plot(train_sizes*100, train_mean, label='average train score (99.99%)', color='blue', marker='o')
    plt.fill_between(train_sizes*100, train_mean + train_std,train_mean - train_std, color='blue', alpha=alpha)
    plt.plot(train_sizes*100, test_mean, label='average validation score (99.88%)', color='red', marker='o')
    plt.fill_between(train_sizes*100, test_mean + test_std, test_mean - test_std, color='red', alpha=alpha)
    plt.xlabel('% of Training Data', fontsize=16)
    plt.ylabel('Accuracy', fontsize=16)
    leg = plt.legend(bbox_to_anchor=(1,0.2), loc="best")
    frame = leg.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    # Plot Annotation
    plt.text(64,0.6,f"Number of Trees: {nTrees}\nTree Depth: {nDepth}\nMin Samples per Leaf: {nLeaf}", size=10, rotation=0.,ha="left", va="center",
         bbox=dict(boxstyle="round",
                   ec='black',
                   fc="white"
                   )
         )
    plt.grid()
    plt.savefig(f"output/learning_curve_dt10_2.jpg", dpi=400)

# Read in the training set
X_train = pd.read_csv('../07_preprocess/model/X_train_dt10_4.csv').iloc[:, :-2]
Y_train = pd.read_csv('../07_preprocess/model/X_train_dt10_4.csv').iloc[:, -2]

# Optimized Parameters that we got from the grid search
nestimators = 10
maxdepth = 3
leaf = 4021

rf_pipeline = Pipeline(steps = [
    ("standardize", StandardScaler()),
    ("rf", RandomForestClassifier(random_state=9, max_depth=maxdepth,  n_estimators=nestimators, min_samples_leaf=leaf))
])

# StratifiedKFold with 10 splits - DO NOT SHUFFLE THE DATA
rkf = StratifiedKFold(n_splits=10, shuffle=False)

# Percentages of data we want to look at
ts = np.linspace(0.1, 1.0, 10)
train_sizes, train_scores, test_scores = learning_curve(estimator=rf_pipeline, X=X_train, y=Y_train.values.ravel(), train_sizes=ts, cv=rkf, scoring='accuracy', n_jobs=-1, verbose=4, shuffle=False)
plot_learning_curve(ts, train_scores, test_scores, nestimators, maxdepth, leaf)
rf_pipeline.fit(X_train, Y_train)

# Get out the final model
joblib.dump(rf_pipeline, 'model/final_model.pkl')

# Save the train and the test scores so you can plot it later
np.savetxt("output/train_scores_finalmodel.txt", train_scores)
np.savetxt("output/test_scores_finalmodel.txt", test_scores)