'''
author: Moshood Olawale
credit: polycode tutorial II
project: dialy_coding drill 01-15-2019 
'''
import numpy as np


class NeuralNetwork():
    def __init__(self):
        np.random.seed(1)
        self.synaptic_weights = 2 * np.random.random((3,1)) - 1


    def sigmoid(self, x):
        return 1/(1+np.exp(-x))

    def sigmoid_derivative(self, x):
        return x*(1-x)

    def train(self, training_inputs, training_outputs, training_iterations):
        for iteration in range(training_iterations):
            outputs = self.think(training_inputs)
            error = training_outputs - outputs
            adjustments = np.dot(training_inputs.T, error*self.sigmoid_derivative(outputs))
            self.synaptic_weights += adjustments  

    def think(self, inputs):
        inputs = inputs.astype(float)
        outputs = self.sigmoid(np.dot(inputs, self.synaptic_weights))
        return outputs


if __name__ == "__main__":
    neural_network = NeuralNetwork()

    print("Random synaptic weights:")
    print(neural_network.synaptic_weights)

    training_inputs = np.array([[0,0,1],
                                [1,1,1],
                                [1,0,1],
                                [0,1,1]])

    training_outputs = np.array([[0,1,1,0]]).T

    neural_network.train(training_inputs, training_outputs, 1000)

    print("Synpatic weights after training:")
    print(neural_network.synaptic_weights)

    A = str(input("Input 1: "))
    B = str(input("Input 2: "))
    C = str(input("Input 3: "))

    print("New testing outputs: ", A, B, C)
    print("Ouput data: ")
    print(neural_network.think(np.array([A,B,C])))