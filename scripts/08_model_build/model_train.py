from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import joblib
import pandas as pd


#Read in the training data
X_train = pd.read_csv('../07_preprocess/model/X_train_dt10_4.csv').iloc[:, :-2]
Y_train = pd.read_csv('../07_preprocess/model/X_train_dt10_4.csv').iloc[:, -2]

# Seeds
seeds = [9, 47, 512]
number_of_splits = 10
number_of_repeats = 3

#
RF_param_grid = {"rf__n_estimators" : [10, 25, 50, 75, 100, 200, 1000], # higher precision, accuracy - reduced overfitting (500,1000) https://bradleyboehmke.github.io/HOML/random-forest.html#
               "rf__max_depth" : [1, 2, 3, 4, 5, 10, None], # https://towardsdatascience.com/mastering-random-forests-a-comprehensive-guide-51307c129cb1#:~:text=max_depth%3A%20The%20number%20of%20splits,shown%20to%20each%20decision%20tree.
             "rf__min_samples_leaf": [1, int(X_train.shape[0] * 0.1/number_of_splits), int(X_train.shape[0] * 0.25/number_of_splits)]
              }

for iteration in range(number_of_repeats):
    print(f"Iteration {iteration+1}")
    rf_pipeline = Pipeline(steps=[
        ("standardize", StandardScaler()),
        ("rf", RandomForestClassifier(random_state=seeds[iteration]))
    ])

    rkf = StratifiedKFold(n_splits=number_of_splits, shuffle=False)
    rfGridSearch = GridSearchCV(rf_pipeline, param_grid=RF_param_grid, scoring='accuracy', refit=False, n_jobs=-1, cv=rkf, verbose=4, return_train_score=True) # 3 iteraitions
    rfGridSearch.fit(X_train, Y_train.values.ravel())
    joblib.dump(rfGridSearch, f'model/gridsearch_final_stratified_{iteration}.pkl')
