from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold, TimeSeriesSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import make_scorer, f1_score

import joblib
import pandas as pd

# apply timeseris split

#Read in the training data
df = pd.read_csv('../07_preprocess/model/X_train_dt10_4.csv')

# make sure all classes have the same samples numbers
minframe= 36320
df = df[df['frameNum'] >= minframe]

X_train = df.iloc[:, :-2]
Y_train = df.iloc[:, -2]

# largest tau * 4
gap = 18160 * 4

# Seeds
seeds = [9, 47, 512]
number_of_splits = 10
number_of_repeats = len(seeds)


#
RF_param_grid = {"rf__n_estimators" : [5, 10, 25, 50, 200, 1000], # higher precision, accuracy - reduced overfitting (500,1000) https://bradleyboehmke.github.io/HOML/random-forest.html#
               "rf__max_depth" : [1, 3 , 5, 10, None], # https://towardsdatascience.com/mastering-random-forests-a-comprehensive-guide-51307c129cb1#:~:text=max_depth%3A%20The%20number%20of%20splits,shown%20to%20each%20decision%20tree.
                 "rf__min_samples_leaf": [1, int((X_train.shape[0]-gap) * 0.1/number_of_splits), int((X_train.shape[0]-gap) * 0.25/number_of_splits)]
              }

# Define scoring dictionary with multiple metrics
scoring = {
    'accuracy': 'accuracy',  # Built-in metric
    'macro_f1': make_scorer(f1_score, average='macro'),  # Custom scorer
}


for iteration in range(number_of_repeats):
    print(f"Iteration {iteration+1}")
    rf_pipeline = Pipeline(steps=[
        ("standardize", StandardScaler()),
        ("rf", RandomForestClassifier(random_state=seeds[iteration]))
    ])


#    rkf = StratifiedKFold(n_splits=number_of_splits, shuffle=False)
    tscv = TimeSeriesSplit(n_splits=number_of_splits,gap=gap,test_size=int((X_train.shape[0]-gap)/(number_of_splits+1))) # test_size needs to be reset bc default doesn't support gap
    rfGridSearch = GridSearchCV(rf_pipeline, param_grid=RF_param_grid, scoring=scoring, refit=False, n_jobs=66, cv=tscv, verbose=4, return_train_score=True) # 3 iteraitions
    rfGridSearch.fit(X_train, Y_train.values.ravel())
    joblib.dump(rfGridSearch, f'model/gridsearch_final_timeseries_t10_{iteration}.pkl')
