import numpy as np
import sympy as sy
from scipy.misc import derivative


def Gauss(A: np.ndarray, b: np.ndarray):
    """
    Implement Gaussian method here to solve A x = b with x unknown.

    Args:
        A: a square matrix defined as a np.array
        b: a vector
    Returns:
        x: the solution vector as a numpy array.
            If the solution does not exists then x should be a
            numpy array with the same dimension as b with np.NaN elements.
            If there is there is infinitly many solutions, x should be the
            solution which set the last free coordinate to 0.
    """
    # DO NOT REMOVE THE FOLLOWING LINES
    # determine if A is multi dimension array
    assert isinstance(A, np.ndarray)
    assert isinstance(b, np.ndarray)
    # determine if len of A is 2 dimension and b is 1 demension
    assert len(A.shape) == 2
    assert len(b.shape) == 1
    # make sure they have same size
    assert b.shape[0] == A.shape[0]
    # TODO
    # Please print: 'No solutions.' if they are no solution
    # Please print: 'Infinitely many solutions.' if they are infinitely many solutions.
    # reference: geeksforgeeks and stackoverflow
    size = len(A[0])
    print(b)
    # elimination
    for row in range(0, size-1):
        for i in range(row+1, size):
            temp = A[i,row] / A[row,row]
            for j in range(row, size):
                A[i,j] = A[i,j] - temp * A[row,j]
            b[i] = b[i] - temp * b[row]

    #print(b)
    #print(A)
    #back tracking
    x = np.zeros((size,1), np.int32)

    # if no solution
    if b[size - 1] != 0 and A[size - 1][size - 1] == 0:
        print("No solutions.")
        for no in x:
            no = np.NaN
        return x

    # if infinite solution
    counter = 0
    for f in range(size-1, 0, -1):
        if b[f] == 0 and A[f][size - 1] == 0:
            counter += 1
        else:
            break;
    if counter != 0:
        for fu in range(0, counter)     :
            x[fu] = 0;

        for row in range(size - 1 - counter, -1, -1):
            sum = b[row]
            for j in range(row + 1, size):
                sum = sum - A[row, j] * x[j]
            x[row] = sum / A[row, row]

        return x

    x[size - 1] = b[size - 1] / A[size - 1, size - 1]
    for row in range(size - 2, -1, -1):
        sum = b[row]
        for j in range(row + 1, size):
            sum = sum - A[row, j] * x[j]
        x[row] = sum / A[row, row]

    return x


def NetwonRaphson(fun, x1: float, x2: float, epsilon=1e-7, x0=0, max_iter=100):
    """
    Netwon Raphson method

    Args:
        fun: the function to find the root (as def f : x --> f(x) real)
        x1, x2: the interval for the search
        epsilon: the update termination condition $f(x) \leq$ epsilon. Default = 1e-7.
        x0: the initial point. Default = 1e-7.
        max_iter: the maximum number of iteration of the Netwon method. Default = 100.
    Returns:
        point: the root of the function if there is a root. None if there is no root.
    """
    # DO NOT REMOVE THE FOLLOWING LINE
    assert hasattr(fun, '__call__'), "fun argument is not a functor!"
    # TODO

    # first get the derivative
    #f_de = fun.diff(x);

    #derivative(f, 5, dx=1e-6)
    # do the iteration
    result = 0.
    for i in range (0, max_iter):
        if result > x1 or result < x2:
            if derivative(fun, result, dx=1e-6) - 0 < epsilon:
                return result
            result = result - fun(result)/derivative(fun, result, dx=1e-6)
        else:
            return result

    if fun(result) > epsilon:
        print("not contain any root")
        exit(0)
    return result


def Netwon(fun, x1: float, x2: float, epsilon=1e-7, x0=0, max_iter=100):
    """
    Netwon Raphson method on the gradient to find optimum point

    Args:
        fun: the function to find the root (as def f : x --> f(x) real)
        x1, x2: the interval for the search
        epsilon: the update termination condition $f(x) \leq$ epsilon. Default = 1e-7.
        x0: the initial point. Default = 1e-7.
        max_iter: the maximum number of iteration of the Netwon method. Default = 100.
    Returns:
        optimum: a point between x1 and x2 that cancel the approximate derivative of fun.
            This point satistfies the first order optimality condition. None if no point cancelled
            the derivative function in the interval.
        type_of_optimum: 1 if it is a local maximum, 0 if it is an inflection point, -1 if it is a local minimum.
            None if the optimum is None.
    """
    # DO NOT REMOVE THE FOLLOWING LINE
    assert hasattr(fun, '__call__'), "fun argument is not a functor!"
    x = NetwonRaphson(fun, x1, x2, epsilon, x0, max_iter)

    x_nega = first_deri(fun, x, -0.01)
    x_posi = first_deri(fun, x, 0.01)
    if x_nega > 0 and x_posi < 0:
        print("maximum")
    if x_nega < 0 and x_posi > 0:
        print("minimum")
    else:
        print("saddle")
    pass

def first_deri(fun, x, small_num):
    return (fun(x + small_num) - fun(x))/small_num
def fun1(x):
    return x * x + 2 * x + 1
def fun2(x):
    return x * x - 1
def fun_test(x):
    return x ** 3 - 5 * x - 9

def main():
    print("question 1")
    A = np.array([[1, 2, 3], [4, 5, 6], [7, 9, 8]], np.int32)
    b = np.array([14, 32, 49], np.int32)
    x = Gauss(A, b)
    print(x)
    print("question 2")
    A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], np.int32)
    b = np.array([14, 32, 49], np.int32)
    x = Gauss(A, b)
    print(x)
    print("question 3")
    A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], np.int32)
    b = np.array([14, 32, 50], np.int32)
    x = Gauss(A, b)
    print(x)

    print("Newton-Paphson test case 1")
    x = NetwonRaphson(fun1, -5, 5, epsilon = 1e-7, x0 = 0, max_iter = 100)
    print(x)
    print("Newton-Paphson test case 2")
    x = NetwonRaphson(fun2, -5, 5, epsilon=1e-7, x0=0, max_iter=100)
    print(x)

    print("Newton test case 1")
    x = Netwon(fun1, -5, 5, epsilon = 1e-7, x0 = 0, max_iter = 100)
    print("Newton test case 2")
    x = Netwon(fun2, -5, 5, epsilon=1e-7, x0=0, max_iter=100)
if __name__ == '__main__':
    main()