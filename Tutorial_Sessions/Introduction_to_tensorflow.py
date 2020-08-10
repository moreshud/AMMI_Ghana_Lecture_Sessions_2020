import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tensorflow as tf

model = tf.keras.Sequential([tf.keras.layers.Dense(units=1, input_shape=[1])])
model.compile("sgd","mean_squared_error")

xs = np.array([-1,0,1,2,3,4], dtype=float)
ys = np.array([-3,-1,1,3,5,7], dtype=float)

model.fit(xs, ys, epochs=500)
print(model.predict(np.array([-29,3,5,10,50.24])))

fashion_mnist = tf.keras.datasets.fashion_mnist.load_data()
(train_images, train_labels), (test_images, test_labels) = fashion_mnist
print(train_images.shape, train_labels.shape)
plt.imshow(train_images[42])
print(train_images[0])
print(train_labels[0])

#normalized the data
train_images, test_images = train_images/255.0, test_images/255.0

# Callback function
class callBack(tf.keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs={}):
        if(logs.get('loss') < 0.4):
            print("\nLoss is low so cancelling training!")
            self.model.stop_training = True

callbacks = callBack()

#building the nn with data
model = tf.keras.models.Sequential([
    tf.keras.layers.Flatten(input_shape=(28, 28)),
    tf.keras.layers.Dense(128, activation=tf.nn.relu),
    tf.keras.layers.Dense(10, activation=tf.nn.softmax)
    ])

# Training the model
model.compile(optimizer=tf.optimizers.Adam(), loss="sparse_categorical_crossentropy")
model.fit(train_images, train_labels, epochs=5, callbacks=[callbacks])
print(model.evaluate(test_images, test_labels))
print(model.predict(test_images))

# Implementing of convolutional and pooling layer
model = tf.keras.models.Sequential([
    tf.keras.layers.Conv2D(64, (3,3), activation="relu", input_shape=(28,28,1)),
    tf.keras.layers.MaxPooling2D(2,2),
    tf.keras.layers.Conv2D(64, (3,3), activation="relu"),
    tf.keras.layers.MaxPooling2D(2,2),    
    tf.keras.layers.Flatten(),
    tf.keras.layers.Dense(128, activation="relu"),
    tf.keras.layers.Dense(10, activation="softmax")
    ])

# summary of the model architecture
model.summary()

# Reshaping the data before training
train_images=train_images.reshape(60000, 28, 28, 1)
train_images=train_images / 255.0
test_images = test_images.reshape(10000, 28, 28, 1)
test_images=test_images/255.0

# Training the model
model.compile(optimizer=tf.optimizers.Adam(), loss="sparse_categorical_crossentropy")
model.fit(train_images, train_labels, epochs=5, callbacks=[callbacks])
print(model.evaluate(test_images, test_labels))
print(model.predict(test_images))