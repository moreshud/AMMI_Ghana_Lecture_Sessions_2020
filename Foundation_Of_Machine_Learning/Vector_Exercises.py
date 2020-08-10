import math

"""
    This exercisises is used to test understanding of Vectors. YOU are NOT to use any Numpy
    implementation for this exercises. 
"""

class Vector(object):
    """
        This class represents a vector of arbitrary size.
        You need to give the vector components. 
    """

    def __init__(self, components=None):
        """
            input: components or nothing
            simple constructor for init the vector
        """
        if components is None:
            components = []
        self.__components = list(components)


    def component(self, i):
        """
            input: index (start at 0)
            output: the i-th component of the vector.
        """
        if type(i) is int and -len(self.__components) <= i < len(self.__components):
            return self.__components[i]
        else:
            raise Exception("index out of range")

    def __len__(self):
        """
            returns the size of the vector
        """
        return len(self.__components)

    def modulus(self):
        """
            returns the euclidean length of the vector
        """
        summe = 0
        
        ## BEGIN SOLUTION
        for i in range(self.__len__()):
            summe += self.__components[i]**2
        return math.sqrt(summe) ## EDIT THIS
        ## END SOLUTION

    def add(self, other):
        """
            input: other vector
            assumes: other vector has the same size
            returns a new vector that represents the sum.
        """
        add_con = []
        size = len(self)
        if size == len(other):
            ## BEGIN SOUTION
            for i in range(size):
                add_con.append(self.__components[i] + other.component(i))
            
            return Vector(add_con) ## EDIT THIS
            ## END SOLUTION
        else:
            raise Exception("must have the same size")

    def sub(self, other):
        """
            input: other vector
            assumes: other vector has the same size
            returns a new vector that represents the difference.
        """
        sub_con = []
        size = len(self)
        if size == len(other):
            ## BEGIN SOUTION
            for i in range(size):
                sub_con.append(self.__components[i] - other.component(i))
            return Vector(sub_con)  ## EDIT THIS
            ## END SOLUTION
        else:  # error case
            raise Exception("must have the same size")

    def multiply(self, other):
        """
            multiply implements the scalar multiplication 
            and the dot-product
        """
        size = len(self)
        if isinstance(other, float) or isinstance(other, int): #scalar multiplicatioj
            ## BEGIN SOLUTION
            sca_con = []
            for i in range(size):
                sca_con.append(other*self.__components[i])            
            return Vector(sca_con) ## EDIT THIS
            ## END SOLUTION
        elif isinstance(other, Vector) and (len(self) == len(other)): # dot product            
            summe = 0
            ## BEGIN SOLUTION
            
            for i in range(size):
                summe += self.__components[i] * other.component(i)
            return summe
            ## END SOLUTION
        else:  # error case
            raise Exception("invalid operand!")

    
    def scalar_proj(self, other):
        """ 
            Computes the scalar projection of vector r on s.
        """
        multiplier = self.multiply(other)
        scalar_con = []
        for i in range(multiplier):
            scalar_con.append(multiplier[i]/self.modulus())
        ### BEGIN SOLUTION
        return scalar_con## EDIT THIS
        ### END SOLUTION
        
    def vector_proj(self, other):
        """ 
            Computes the vector projection of vector r on s.
            use the other functions created above.
        """
    
        ### BEGIN SOLUTION
        multiplier = self.multiply(other)
        vector_con = []
        for i in range(other.__len__()):
            vector_con.append((multiplier/self.modulus()**2)*other.component(i))
        
        return vector_con ## EDIT THIS
        ### END SOLUTION


##################
## TESTING CLASS ######
##################
     
import unittest

class Test(unittest.TestCase):
    def test_component(self):
        """
            test for method component
        """
        x = Vector([1, 2, 3])
        self.assertEqual(x.component(0), 1)
        self.assertEqual(x.component(2), 3)
        try:
            y = Vector()
            self.assertTrue(False)
        except:
            self.assertTrue(True)


    def test_size(self):
        """
            test for size()-method
        """
        x = Vector([1, 2, 3, 4])
        self.assertEqual(len(x), 4)

    def test_modulus(self):
        """
            test for the eulidean length
        """
        x = Vector([1, 2])
        self.assertAlmostEqual(x.modulus(), 2.236, 3)

    def test_add(self):
        """
            test for + operator
        """
        x = Vector([1, 2, 3])
        y = Vector([1, 1, 1])
        self.assertEqual((x.add(y)).component(0), 2)
        self.assertEqual((x.add(y)).component(1), 3)
        self.assertEqual((x.add(y)).component(2), 4)

    def test_sub(self):
        """
            test for - operator
        """
        x = Vector([1, 2, 3])
        y = Vector([1, 1, 1])
        self.assertEqual((x.sub(y)).component(0), 0)
        self.assertEqual((x.sub(y)).component(1), 1)
        self.assertEqual((x.sub(y)).component(2), 2)

    def test_multiply(self):
        """
            test for vector multiplication
        """
        x = Vector([1, 2, 3])
        a = Vector([2, -1, 4])  # for test of dot-product
        b = Vector([1, -2, -1])
        self.assertEqual((x.multiply(3.0)).component(0),3 )
        self.assertEqual((a.multiply(b)), 0)

    def test_scalar_projection(self):
        """
            test for scalar projection
        """
        x = Vector([3, 4])
        y = Vector([4, 3])
        self.assertEqual(x.scalar_proj(y), 4.8)
        
    def test_vector_projection(self):
        """
            test for scalar projection
        """
        x = Vector([3, 4])
        y = Vector([4, 3])
        self.assertEqual((x.vector_proj(y)).component(1), 3.84)
        
        

if __name__ == "__main__":
    unittest.main()

