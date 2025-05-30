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

#def time_series_learning_curve(estimator, X, y, cv_splits=20, train_sizes=np.linspace(0.1, 1.0, 10), scoring=accuracy_score):
#    """
#    Custom implementation of learning curve for TimeSeriesSplit.
#
#    Parameters:
#    - estimator: Sklearn model pipeline (must implement fit & predict).
#    - X: Feature matrix.
#    - y: Target vector.
#    - cv_splits: Number of splits for TimeSeriesSplit.
#    - train_sizes: List of proportions of data to use for training.
#    - scoring: Scoring function (default: accuracy_score).
#
#    Returns:
#    - train_sizes_abs: Actual training sizes.
#    - train_scores: Array of shape (len(train_sizes), cv_splits) with training scores per fold.
#    - test_scores: Array of shape (len(train_sizes), cv_splits) with validation scores per fold.
#    """
#    tscv = TimeSeriesSplit(n_splits=cv_splits)
#    max_train_size = len(X) - (len(X) // (cv_splits + 1))  # Compute max train size
#    train_sizes_abs = (train_sizes * max_train_size).astype(int)  # Ensure valid training sizes
#
#    train_scores = []
#    test_scores = []
#
#
#    
#    for train_size in train_sizes_abs:
#        fold_train_scores = []
#        fold_test_scores = []
#
#        for train_index, test_index in tscv.split(X):
#            if len(train_index) < train_size:
#                print(f"training data {len(train_index)} smaller than train_size {train_size}")
#                X_train, X_test = X.iloc[train_index], X.iloc[test_index]
#                y_train, y_test = y.iloc[train_index], y.iloc[test_index]
#                
#                
#            X_train, X_test = X.iloc[train_index[:train_size]], X.iloc[test_index]
#            y_train, y_test = y.iloc[train_index[:train_size]], y.iloc[test_index]
#
#            estimator.fit(X_train, y_train)
#            y_train_pred = estimator.predict(X_train)
#            y_test_pred = estimator.predict(X_test)
#
#            fold_train_scores.append(scoring(y_train, y_train_pred))
#            fold_test_scores.append(scoring(y_test, y_test_pred))
#
#        train_scores.append(fold_train_scores)  # Store all fold scores
#        test_scores.append(fold_test_scores)
#
#    return train_sizes_abs, np.array(train_scores), np.array(test_scores)

def time_series_learning_curve(estimator, X, y, cv_splits=20,gap=10, train_sizes=np.linspace(0.1, 1.0, 10), scoring=accuracy_score):
    """
    Custom implementation of learning curve for TimeSeriesSplit.

    Parameters:
    - estimator: Sklearn model pipeline (must implement fit & predict).
    - X: Feature matrix.
    - y: Target vector.
    - cv_splits: Number of splits for TimeSeriesSplit.
    - train_sizes: List of proportions of data to use for training.
    - scoring: Scoring function (default: accuracy_score).

    Returns:
    - train_sizes_abs: Actual training sizes.
    - train_scores: Array of shape (len(train_sizes), cv_splits) with training scores per fold.
    - test_scores: Array of shape (len(train_sizes), cv_splits) with validation scores per fold.
    """
    tscv = TimeSeriesSplit(n_splits=cv_splits,gap=gap,test_size=int((X.shape[0]-gap)/(cv_splits+1)))
    max_train_size = len(X) - (len(X) // (cv_splits + 1))  # Compute max train size
    #train_sizes_abs = (train_sizes * max_train_size).astype(int)  # Ensure valid training sizes
    
    train_sizes_abs=[]
    train_scores = []
    test_scores = []
    

    

    for train_index, test_index in tscv.split(X):
                
        train_sizes_abs.append(len(train_index))
        
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        estimator.fit(X_train, y_train)
        y_train_pred = estimator.predict(X_train)
        y_test_pred = estimator.predict(X_test)
        
        train_scores.append(scoring(y_train, y_train_pred))
        test_scores.append(scoring(y_test, y_test_pred))

    return np.array(train_sizes_abs), np.array(train_scores)[:,None], np.array(test_scores)[:,None]


def plot_learning_curve(train_sizes,frame_tot, train_scores, test_scores, params, save_path="output/learning_curve.jpg",alpha=0.2):
    plt.figure(figsize=(8, 6))

    # Compute mean & standard deviation
    train_mean, train_std = np.mean(train_scores, axis=1), np.std(train_scores, axis=1)
    test_mean, test_std = np.mean(test_scores, axis=1), np.std(test_scores, axis=1)

    # Plot Training & Validation Scores
    plt.plot(train_sizes /frame_tot * 100, train_mean, 'o-', label="Training Accuracy", color='blue')
    plt.fill_between(train_sizes /frame_tot * 100, train_mean - train_std, train_mean + train_std, color='blue', alpha=alpha)

    plt.plot(train_sizes/frame_tot * 100, test_mean, 'o-', label="Validation Accuracy", color='red')
    plt.fill_between(train_sizes/frame_tot * 100, test_mean - test_std, test_mean + test_std, color='red', alpha=alpha)

    # Labels and Annotations
    plt.xlabel('% of Training Data', fontsize=14)
    plt.ylabel('Accuracy', fontsize=14)
    plt.legend(loc="best")
    plt.title('Learning Curve')
    plt.grid()

    # Get axis limits dynamically
    x_min, x_max = plt.xlim()
    y_min, y_max = plt.ylim()
    # Display Model Parameters
    param_text = f"Number of Trees: {params['n_estimators']}\nDepth: {params['max_depth']}\nMin Leaf Samples: {params['min_samples_leaf']}"
    plt.text(0.65*x_max, 0.6*y_max, param_text, bbox=dict(boxstyle="round", fc="white", ec="black"))

#    # Plot Annotation
#    plt.text(64,0.6,f"Number of Trees: {nTrees}\nTree Depth: {nDepth}\nMin Samples per Leaf: {nLeaf}", size=10, rotation=0.,ha="left", va="center",
#         bbox=dict(boxstyle="round",
#                   ec='black',
#                   fc="white"
#                   )
#         )


    # Save Plot
    plt.savefig(save_path, dpi=400)
    plt.show()


# Read in the training set
df = pd.read_csv('../07_preprocess/model/X_train_dt10_4.csv')
X_train = df.iloc[:, :-2]
Y_train = df.iloc[:, -2]

sample_tot = df.shape[0]
print(f"total frame number: {sample_tot}")

# largest tau * 4
gap = 18160 * 4

# Optimized Parameters that we got from the grid search
nestimators = 200
maxdepth = 3
leaf = 1641

#rf_pipeline = Pipeline(steps = [
#    ("standardize", StandardScaler()),
#    ("rf", RandomForestClassifier(random_state=9, max_depth=maxdepth,  n_estimators=nestimators, min_samples_leaf=leaf))
#])

rf_pipeline = Pipeline([
    ("scaler", StandardScaler()), 
    ("rf", RandomForestClassifier(
        random_state=9, 
        max_depth=maxdepth,  
        n_estimators=nestimators, 
        min_samples_leaf=leaf, 
        n_jobs=30,  # Added for parallel processing
    ))
])



train_sizes = np.linspace(0.1, 1.0, 10)
# === Manually Set Larger Training Sizes ===
#train_sizes = np.linspace(0.01, 1.0, 10) * max_train_size  # 1% to 100% of dataset
#train_sizes = train_sizes.astype(int)  # Convert to integer values

# learning curve calculation
#train_sizes, train_scores, test_scores = learning_curve(estimator=rf_pipeline, X=X_train, y=Y_train.values.ravel(), train_sizes=train_sizes, cv=rkf, scoring='accuracy', n_jobs=36, verbose=4, shuffle=False)

cv_splits=20

# === Compute Custom Learning Curve ===
train_sizes_abs, train_scores, test_scores = time_series_learning_curve(
    estimator=rf_pipeline, 
    X=X_train, 
    y=Y_train, 
    cv_splits=cv_splits,
    gap=gap
)



# plot learning curve
plot_learning_curve(train_sizes_abs,(sample_tot-gap)*cv_splits/(cv_splits+1), train_scores, test_scores, {"n_estimators": nestimators, "max_depth": maxdepth, "min_samples_leaf": leaf})

#plot_learning_curve(train_sizes, train_scores, test_scores, nestimators, maxdepth, leaf)

## train the final model
rf_pipeline.fit(X_train, Y_train)

# Get out the final model
joblib.dump(rf_pipeline, 'model/final_model_0210.pkl')

# Save the train and the test scores so you can plot it later
np.savetxt("output/train_scores_finalmodel_timeseries.txt", train_scores)
np.savetxt("output/test_scores_finalmodel_timeseries.txt", test_scores)
