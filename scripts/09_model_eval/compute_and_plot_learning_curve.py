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

    # Display Model Parameters
    param_text = f"Number of Trees: {params['n_estimators']}\nDepth: {params['max_depth']}\nMin Leaf Samples: {params['min_samples_leaf']}"
    plt.text(65, 0.6, param_text, bbox=dict(boxstyle="round", fc="white", ec="black"))

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



# StratifiedKFold with 10 splits - DO NOT SHUFFLE THE DATA
rkf = StratifiedKFold(n_splits=20, shuffle=False)
#tscv = TimeSeriesSplit(n_splits=5,gap=gap,test_size=int((X_train.shape[0]-gap)/(number_of_splits+1))) # test_size needs to be reset bc default doesn't support gap
#tscv = TimeSeriesSplit(n_splits=10) # a simpler version of timeseries split
# === Get Maximum Training Size Allowed ===
#max_train_size = sample_tot - (sample_tot // (tscv.n_splits + 1))

#cv_iter = tscv.split(X_train)

# Percentages of data we want to look at
train_sizes = np.linspace(0.1, 1.0, 20)
# === Manually Set Larger Training Sizes ===
#train_sizes = np.linspace(0.01, 1.0, 10) * max_train_size  # 1% to 100% of dataset
#train_sizes = train_sizes.astype(int)  # Convert to integer values

# learning curve calculation
train_sizes, train_scores, test_scores = learning_curve(estimator=rf_pipeline, X=X_train, y=Y_train.values.ravel(), train_sizes=train_sizes, cv=rkf, scoring='accuracy', n_jobs=36, verbose=4, shuffle=False)

print(f"train_sizes: {train_sizes}")

# plot learning curve
plot_learning_curve(train_sizes,sample_tot, train_scores, test_scores, {"n_estimators": nestimators, "max_depth": maxdepth, "min_samples_leaf": leaf})

#plot_learning_curve(train_sizes, train_scores, test_scores, nestimators, maxdepth, leaf)

## train the final model
rf_pipeline.fit(X_train, Y_train)

# Get out the final model
joblib.dump(rf_pipeline, 'model/final_model_0210.pkl')

# Save the train and the test scores so you can plot it later
np.savetxt("output/train_scores_finalmodel.txt", train_scores)
np.savetxt("output/test_scores_finalmodel.txt", test_scores)
