
import tensorflow as tf

from keras import optimizers
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation, Flatten, Reshape, Permute,Conv1D,MaxPooling1D,Input,Subtract,Multiply,Add,Lambda,concatenate
from keras.callbacks import LearningRateScheduler
import keras.backend as K
from keras.callbacks import EarlyStopping
from keras.models import Model


l1=84
k=2
inputX = Input(shape=(l1,1))
inputPave = Input(shape=(l1,1)) #pave
inputNave = Input(shape=(l1,1))

x = Conv1D(16, kernel_size=k, activation='relu')(inputX)
x = Conv1D(32, kernel_size=k, activation='relu')(x)
x = Flatten()(x)

p = Conv1D(16, kernel_size=k, activation='relu')(inputPave)
p = Conv1D(32, kernel_size=k, activation='relu')(p)
p = Flatten()(p)

n = Conv1D(16, kernel_size=k, activation='relu')(inputNave)
n = Conv1D(32, kernel_size=k, activation='relu')(n)
n = Flatten()(n)

pbia = Subtract()([x,p]) # p.output
nbia = Subtract()([x,n])

dropoutp1 = Dropout(0.2)(pbia)
dropoutp2 = Dropout(0.5)(pbia)
hp1 = Dense(128, activation='relu')(dropoutp1) #
hp2 = Dense(1, activation='sigmoid')(dropoutp2) #
hp3 = Dense(1, activation='sigmoid')(hp1) #
outputp = Add()([hp2, hp3])

dropoutn1 = Dropout(0.2)(nbia)
dropoutn2 = Dropout(0.5)(nbia)
hn1 = Dense(128, activation='relu')(dropoutn1) #
hn2 = Dense(1, activation='sigmoid')(dropoutn2) #
hn3 = Dense(1, activation='sigmoid')(hn1) #

outputn = Add()([hn2, hn3])

output = concatenate([outputn, outputp],axis=1) #

model = Model(inputs=[inputX,inputPave,inputNave], outputs=output)

model.compile(optimizer='adamax',
              loss='mean_squared_error',
              metrics=['accuracy'])
