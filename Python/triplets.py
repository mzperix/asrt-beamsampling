import numpy as np
from ideal_observer import PRED_PROB_ALL
import pandas as pd

HIGH = "high"
LOW = "low"
UNIFORM = "uniform"

class TripletModel():
    def __init__(self,**kwargs):
        self.kwargs = kwargs

    def predict_rt(self,data,data_train=None):
        return(triplet_model_predict(data,data_train,predict_all=False))

    def generate_predictive_probabilities(self,data):
        return({PRED_PROB_ALL:triplet_model_predict(data, predict_all=True)})

def mark_triplet_type(data, sequence, predict_all=False):
    Y = data['Y']
    block = data['block']

    y1 = Y[0]
    y2 = Y[1]
    b1 = block[0]
    b2 = block[1]

    if predict_all:
        predictions = [[UNIFORM]*4,
                       [UNIFORM]*4]
    else:
        predictions = [UNIFORM, UNIFORM]
    for j,y in enumerate(Y[2:]):
        i = np.where(sequence == y1)[0][0]
        if block[j] == b1:
            if predict_all:
                prediction = [UNIFORM, UNIFORM, UNIFORM, UNIFORM]
                prediction[sequence[i+1]] = HIGH
                predictions.append(prediction)
            else:
                if sequence[i+1]==y:
                    predictions.append(HIGH)
                else:
                    predictions.append(LOW)

        else:
            if predict_all:
                predictions.append([UNIFORM]*4)
            else:
                predictions.append(UNIFORM)

        y1 = y2
        y2 = y
        b1 = b2
        b2 = block[j]
    return predictions

def triplet_model_predict(data, data_train=None, predict_all=False):
    if data_train is None:
        data_train = data
    sequence = np.array(data_train['Y'][data_train['trial_type']=='P'][:5])
    prediction_dict = {
        HIGH: np.log(0.625),
        LOW: np.log(0.125),
        UNIFORM: np.log(0.25)
        }
    tags = mark_triplet_type(
        data=data,
        sequence=sequence,
        predict_all=predict_all)

    if predict_all:
        predictions = [[prediction_dict[p] for p in row] for row in tags]
        return(np.array([np.transpose(np.array(predictions))]))
    else:
        tags_train = mark_triplet_type(
            data=data_train,
            sequence=sequence,
            predict_all=predict_all
        )
        df = pd.DataFrame(
            dict(
                rt=data_train['rt'],
                triplet_type=tags_train,
                correct_response=data_train['correct_response'],
                filters=data_train['filters'])
        )
        means = {}
        for k, v in prediction_dict.items():
            means[k] = df[
                (df['correct_response']==1) &
                (df['filters']==True) &
                (df['triplet_type']==k)]\
                .mean()['rt']
        predictions = [means[p] for p in tags]
        return(np.array(predictions))

def triplet_correlation(data):
    return(np.corrcoef(
        triplet_model_predict(data)[data['correct_response']==1],
        data['rt'][data['correct_response']==1])[0,1])

def triplet_overlap(data1,data2):
    sequence1 = np.array(data1['Y'][data1['trial_type']=='P'][:5])
    sequence2 = np.array(data2['Y'][data2['trial_type']=='P'][:5])
    print(sequence1, sequence2)
    counter = 0
    for j,y in enumerate(sequence1[:4]):
        i = np.where(sequence2 == y)[0][0]
        if sequence2[i+1]==sequence1[j+1]:
            counter += 1
    return counter
