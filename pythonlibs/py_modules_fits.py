import sys 
import numpy as np 

"""
Calculates the value of a Gaussian distribution at a given point.

Parameters:
    x (float or numpy.ndarray): The input value(s) at which the Gaussian distribution is evaluated.
    mean (float): The mean of the Gaussian distribution.
    sigma (float): The standard deviation of the Gaussian distribution.

Returns:
    float or numpy.ndarray: The value(s) of the Gaussian distribution evaluated at the given input(s).

"""
def gaussian(x, mean, sigma):
    return (1.0 / (np.sqrt(2.0 * np.pi*sigma**2))) * np.exp(-0.5 * ((x - mean)/sigma)**2)


"""
Calculates the sum of three Gaussians at a given point.

Parameters:
    x (float or numpy.ndarray): The input value(s) at which the sum of Gaussians is evaluated.
    mean1 (float): The mean of the first Gaussian.
    sigma1 (float): The standard deviation of the first Gaussian.
    weight1 (float): The weight for the first Gaussian.
    mean2 (float): The mean of the second Gaussian.
    sigma2 (float): The standard deviation of the second Gaussian.
    weight2 (float): The weight for the second Gaussian.
    mean3 (float): The mean of the third Gaussian.
    sigma3 (float): The standard deviation of the third Gaussian.
    weight3 (float): The weight for the third Gaussian.

Returns:
    float or numpy.ndarray: The value(s) of the sum of three Gaussians evaluated at the given input(s).

"""
def gaussian_2(x, mean1, sigma1, weight1, mean2, sigma2, weight2):

    gaussian_sum_value = (weight1 * gaussian(x, mean1, sigma1) +
                          weight2 * gaussian(x, mean2, sigma2))

    return gaussian_sum_value


"""
Calculates the sum of multiple Gaussian functions at the given x-values.

Parameters:
    x (array-like): The x-values at which the Gaussian functions are evaluated.
    *params: Variable-length argument list containing the parameters for the Gaussian functions.
                The parameters are organized in groups of three: amplitude, mean, and standard deviation.
                The number of parameters should be a multiple of 3.

Returns:
    float: The sum of the Gaussian functions evaluated at the given x-values.
    
"""
def gaussian_X(x,*params):

    gaussian_sum_value = 0
    for i in range(0,int(len(params)),3):
        gaussian_sum_value += params[2+i] * gaussian(x, params[i], params[1+i])

    return gaussian_sum_value


"""
Calculates the value of a linear fit at the given x-values.

Parameters:
    x (float or numpy.ndarray): The input value(s) at which the linear fit is evaluated.
    a (float): Slope of the linear curve.
    b (float): y-direction shift.

Returns:
    float or numpy.ndarray: The value(s) of the linear fit evaluated at the given input(s).

"""
def linear(x,a,b):
    return a*x+b

