from sklearn.model_selection import train_test_split
from tensorflow.keras.optimizers import Adam
import pandas as pd
from tensorflow.keras.layers import Input, Bidirectional, LSTM, Concatenate, Dense, Dropout
from tensorflow.keras.models import Model
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import roc_curve, precision_recall_curve


def model_2(RNA_GE_feature_length, lstm_layer_units=[16, 16], hidden_layer_units=[64, 32], dropout=0.3):

    # Inputs
    input_RNA_left = Input(shape=(RNA_GE_feature_length, 1), dtype='float64')
    input_RNA_right = Input(shape=(RNA_GE_feature_length, 1), dtype='float64')
   
    x_left = input_RNA_left
    x_right = input_RNA_right
   
    # LSTM layers
    for i in range(len(lstm_layer_units) - 1):  # First few LSTM layers
        lstm = Bidirectional(LSTM(lstm_layer_units[i], return_sequences=True))
        x_left = lstm(x_left)
        x_right = lstm(x_right)
   
    # Last LSTM layer
    lstm = Bidirectional(LSTM(lstm_layer_units[-1], return_sequences=False))
    x_left = lstm(x_left)
    x_right = lstm(x_right)
   
    # Concatenate
    x = Concatenate()([x_left, x_right])
    
    # Dense layers for classification
    for i in range(len(hidden_layer_units)):
        x = Dense(units=hidden_layer_units[i], activation="relu")(x)
        if dropout:
            x = Dropout(rate=dropout)(x)
   
    # Output layer
    output = Dense(1, activation='sigmoid')(x)
    
    # Create and return the model
    model = Model(inputs=[input_RNA_left, input_RNA_right], outputs=output)
    return model
    
import os,sys

dirL = os.getcwd() # the path of current script
eachcell = sys.argv[1] # gm12878_ctcf
dirX = sys.argv[2] # out_dir_test
# pathFile = dirL+'/data'# data/ path


filename = f'{eachcell}_distance_matched.csv_winGEfea.csv'

file_path = dirL+f'/{dirX}/{filename}'
print(file_path)

data = pd.read_csv(file_path)
left_features = data.filter(regex='^L').values.reshape(-1, data.filter(regex='^L').shape[1], 1)
right_features = data.filter(regex='^R').values.reshape(-1, data.filter(regex='^R').shape[1], 1)
labels = data['label'].values


# split-data
left_train, left_test, right_train, right_test, y_train, y_test = train_test_split(
    left_features, right_features, labels, test_size=0.2, random_state=42)

model = model_2(RNA_GE_feature_length=left_features.shape[1], lstm_layer_units=[16, 16], hidden_layer_units=[64, 32], dropout=0.3)

# Compile

model.compile(optimizer=Adam(learning_rate=0.001), loss='binary_crossentropy', metrics=['accuracy'])

# tain
epoch = 30
batch_size = 30
history = model.fit([left_train, right_train], y_train, validation_data=([left_test, right_test], y_test), epochs=epoch, batch_size=batch_size)

#model.save(f'{filename}_model.h5') 

loss, accuracy = model.evaluate([left_test, right_test], y_test)
print(f'Test loss: {loss}, Test accuracy: {accuracy}')

y_pred_proba = model.predict([left_test, right_test])
auc = roc_auc_score(y_test, y_pred_proba)
auprc = average_precision_score(y_test, y_pred_proba)
print(f"AUC:{auc}, AUPRC:{auprc}")

fpr, tpr, _ = roc_curve(y_test, y_pred_proba)

precision, recall, _ = precision_recall_curve(y_test, y_pred_proba)

roc_data = pd.DataFrame({'FPR': fpr, 'TPR': tpr})
pr_data = pd.DataFrame({'Precision': precision, 'Recall': recall})


roc_data.to_csv(f'{file_path}.roc_data.csv', index=False)
pr_data.to_csv(f'{file_path}.pr_data.csv', index=False)
