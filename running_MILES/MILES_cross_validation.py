import numpy as np
import tensorflow as tf
from mil.data.datasets.loader import load_data
from mil.data.datasets import musk1, protein
from mil.data.datasets import musk1, protein
from mil.trainer import Trainer
from mil.bag_representation.mapping import DiscriminativeMapping
from mil.preprocessing import StandarizerBagsList
from mil.validators import KFold
from mil.models import MILES
from mil.metrics import BinaryAccuracy, Precision, Recall, AUC
from logzero import logger
from sklearn.utils import resample
import pandas as pd

# v02: 
# no normalisation
# use from mil.data.datasets.loader import load_data

class Test:
    def __init__(self, sigma2, c, in_file = ""):
        self.sigma2 = sigma2
        self.infile = in_file
        self.c = c
    def load(self):
        return(load_data(self.infile))
    def load_dataset(self):
        """Save thing to a file."""
        np.random.seed(42)
        #(self.bags_train, self.y_train),(self.bags_test, self.y_test) = musk1.load()
        (self.bags_train, self.y_train),(self.bags_test, self.y_test) = self.load()
        # full dataset
        self.bags_all = self.bags_train
        self.bags_all.extend(self.bags_test)
        self.y_all = np.concatenate((self.y_train, self.y_test))

    def run_model(self):
        
        self.trainer = Trainer()
        # metrics come from tenserflow library tf.keras.metrics
        metrics = [AUC, BinaryAccuracy, Precision, Recall]
        model = MILES(sigma2=self.sigma2, C=self.c)
        logger.info(self.sigma2)
        self.trainer.prepare(model, metrics=metrics)
        # inner 10fold
        valid = KFold(n_splits=10, shuffle=True)
        seeds = [42,135,23,7,49,1991,1970,30,108,5]
        self.performance = {"auc": [], "binaryaccuracy": [], "precision" : [], "recall": [], "seed": []}
        # outer 10 fold 
        for seed in seeds:
            np.random.seed(seed)
            history = self.trainer.fit(self.bags_all, self.y_all, validation_strategy=valid, verbose=1)
            auc =[e["auc"] for e in history['metrics_val']]
            binaryaccuracy = [e["binaryaccuracy"] for e in history['metrics_val']]
            precision = [e["precision"] for e in history['metrics_val']]
            recall = [e["recall"] for e in history['metrics_val']]
            seed = [seed]*10
            self.performance["auc"].extend(auc)
            self.performance["binaryaccuracy"].extend(binaryaccuracy)
            self.performance["precision"].extend(precision)
            self.performance["recall"].extend(recall)
            self.performance["seed"].extend(seed)
        self.performance = pd.DataFrame(self.performance)





if __name__ == "__main__":

    import dill
    import sys

    sigma2=float(sys.argv[1])
    c=float(sys.argv[2])
    infile=sys.argv[3]
    #model_file = sys.argv[4]
    performance_file = sys.argv[4]
    # code for standalone use
    test = Test(sigma2, c, in_file=infile)
    test.load_dataset()
    test.run_model()

    #dill.dump(test.trainer, open(model_file, 'wb'))
    #obj = dill.load(open(model_file, 'rb'))
    #print(obj.metrics_train)
    test.performance.to_csv(performance_file, sep="\t")
