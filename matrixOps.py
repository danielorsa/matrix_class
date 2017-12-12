from matrixClass import *
import math

# This module employs my matrix class to perform matrix operations

#### Covariance ####
def covariance(matrix0):
    ''' Compute and return a covariance matrix.
    :param matrix0: A Matrix object.
    :return covMatrix: A Matrix object.
    '''
    covMatrix = Matrix()
    # instantiate covariance matrix
    for i in range(matrix0.columns):
        row = []
        for j in range(matrix0.columns):
            row.append(0)
        covMatrix.rows += 1
        covMatrix.columns +=1
        covMatrix.addRow(row)
    k = matrix0.rows
    matrixList = []
    meanMatrix = mean(matrix0)

    for i in range(matrix0.rows):
        # create vector for every data pair
        vector = Matrix()
        vector.columns += 1
        for j in range(matrix0.columns):
            vector.addRow([matrix0.matrix[i][j]])
            vector.rows += 1
        vector = matrixSub(vector, meanMatrix)  # subtract mean vector from measurement vector
        transVector = transpose(vector)
        prodMatrix = matrixMult(vector, transVector) # multiply vector by its transpose
        matrixList.append(prodMatrix)

    for i in matrixList:
        covMatrix = matrixAdd(covMatrix, i) # sum all vector*transpose products

    covMatrix = scalarMult((1/k), covMatrix) # multiply sum vector by 1/k
    return covMatrix


### Transpose ###
def transpose(matrix0):
    ''' Creates a transposed version of a matrix.
    :param matrix0: A Matrix object.
    :return tMatrix: A Matrix object.
    '''
    tMatrix = Matrix()
    tMatrix.rows, tMatrix.columns = matrix0.columns, matrix0.rows
    tMatrix.totalElements = tMatrix.rows * tMatrix.columns
    for i in range(matrix0.columns):
        tMatrix.matrix.append([])
        for j in range(matrix0.rows):
            tMatrix.matrix[i].append(matrix0.matrix[j][i])
    return tMatrix


### Mean ###
def mean(matrix0):
    ''' Computes the mean vector of a matrix.
    :param matrix0: A Matrix object.
    :return meanMatrix: A Matrix object.
    '''
    meanMatrix = Matrix()
    meanMatrix.rows = matrix0.columns
    meanMatrix.columns = 1
    meanMatrix.totalElements = meanMatrix.rows

    for i in range(matrix0.columns):
        row = [0]
        for j in range(matrix0.rows):
            row[0] += matrix0.matrix[j][i]
        meanMatrix.addRow(row)
    meanMatrix = scalarMult((1/matrix0.rows), meanMatrix)
    return meanMatrix


#### Add ####
def matrixAdd(matrix0, matrix1):
    ''' Computes the sum of two conformable matrices.
    :param matrix0: A Matrix object.
    :param matrix1: A Matrix object.
    :return sumMatrix: A Matrix object.
    '''
    sumMatrix = Matrix()
    sumMatrix.rows = matrix0.rows
    sumMatrix.columns = matrix0.columns
    sumMatrix.totalElements = sumMatrix.rows * sumMatrix.columns

    if (matrix0.rows == matrix1.rows) and (matrix0.columns == matrix1.columns):
        for i in range(matrix0.rows):
            row = []
            for j in range(matrix0.columns):
                row.append(matrix0.matrix[i][j] + matrix1.matrix[i][j])
            sumMatrix.addRow(row)
        return sumMatrix


#### Sub ####
def matrixSub(matrix0, matrix1):
    ''' Computes the difference of two conformable matrices.
    :param matrix0: A Matrix object.
    :param matrix1: A Matrix object.
    :return diffMatrix: A Matrix object.
    '''
    diffMatrix = Matrix()
    diffMatrix.rows = matrix0.rows
    diffMatrix.columns = matrix0.columns
    diffMatrix.totalElements = diffMatrix.rows * diffMatrix.columns

    if (matrix0.rows == matrix1.rows) and (matrix0.columns == matrix1.columns):
        for i in range(matrix0.rows):
            row = []
            for j in range(matrix0.columns):
                row.append(matrix0.matrix[i][j] - matrix1.matrix[i][j])
            diffMatrix.addRow(row)
        return diffMatrix


#### Scalar Mult ####
def scalarMult(scalar, matrix0):
    ''' Performs a scalar multiplication on a matrix.
    :param scalar: a numerical data type; int/float.
    :param matrix0: A Matrix object.
    :return sProdMatrix: A Matrix object.
    '''
    sProdMatrix = Matrix()
    sProdMatrix.rows = matrix0.rows
    sProdMatrix.columns = matrix0.columns
    sProdMatrix.totalElements = sProdMatrix.rows * sProdMatrix.columns

    for line in matrix0.matrix:
        row = []
        for element in line:
            row.append(element*scalar)
        sProdMatrix.addRow(row)
    return sProdMatrix


#### Mult ####
def matrixMult(matrix0, matrix1):
    ''' Computes the product of two conformable matrices.
    :param matrix0: A Matrix object.
    :param matrix1: A Matrix object.
    :return prodMatrix: A Matrix object.
    '''
    prodMatrix = Matrix()
    prodMatrix.rows = matrix0.rows
    prodMatrix.columns = matrix1.columns
    prodMatrix.totalElements = prodMatrix.rows * prodMatrix.columns

    if (matrix0.columns == matrix1.rows):

        for i in range(matrix0.rows):
            row = []
            for j in range(matrix1.columns):
                row.append(0)
            prodMatrix.addRow(row)

        for i in range(matrix0.rows):
            for j in range(matrix1.columns):
                for k in range(matrix1.rows):
                    prodMatrix.matrix[i][j] += matrix0.matrix[i][k] * matrix1.matrix[k][j]
        return prodMatrix

    elif (matrix0.rows == matrix1.columns):
        return matrixMult(matrix1, matrix0)


#### Discriminant function ####
def discriminant(vector, dataSet):
    ''' Computes the discriminant of a measurement vector against a data set.
    :param vector: A Matrix object; measurement vector.
    :param dataSet: A Matrix object; set of data.
    :return discrim: A numerical data type; int/float.
    '''
    a = (-.5)*(matrixMult(matrixMult((transpose(matrixSub(vector, mean(dataSet)))), (covariance(dataSet)).inverse()), matrixSub(vector, mean(dataSet)))).matrix[0][0]
    b = (-.5)*math.log((covariance(dataSet)).determinant())
    c = math.log(.5)
    discrim = a + b + c
    return discrim


#### G.J. Elimination ####
def gaussJordan(matrix0):
    ''' Computes the solutions of a system by the Gauss-Jordan method.
    :param matrix0: A Matrix object; augmented linear system.
    :return gjMatrix: A matrix object; solutions to the linear system.
    '''
    gjMatrix = Matrix()
    for i in matrix0.matrix:
        row = []
        for j in i:
            row.append(j)
        gjMatrix.addRow(row)
    gjMatrix.columns = matrix0.columns
    gjMatrix.rows = matrix0.rows

    for j in range(gjMatrix.rows):
        p = 0
        greatest = 0
        for i in range(j, gjMatrix.rows):   # compute pivot index
            if abs(gjMatrix.matrix[i][j]) > greatest:
                greatest = gjMatrix.matrix[i][j]
                p = i
        if gjMatrix.matrix[p][j] == 0:
            break
        if p > j:   # interchange row p and j
            gjMatrix.matrix[p], gjMatrix.matrix[j] = gjMatrix.matrix[j], gjMatrix.matrix[p]
        for i in range(gjMatrix.columns):   # divide row j by C[j][j]
            newElement = gjMatrix.matrix[j][i]/greatest
            gjMatrix.matrix[j][i] = newElement
        for i in range(gjMatrix.rows):  # for each i =/= j, subtract C[i][j] times row j from row i
            eliminator = gjMatrix.matrix[i][j]
            if i != j:
                for x in range(gjMatrix.columns):
                    gjMatrix.matrix[i][x] -= eliminator * gjMatrix.matrix[j][x]

    for i in range(gjMatrix.rows): # partition solutions
        newRow = []
        newRow.append(gjMatrix.matrix[i][-1])
        gjMatrix.matrix[i] = newRow
    gjMatrix.columns = 1
    return gjMatrix


#### Gaussian Elimination ####
def gaussElim(matrix0):
    ''' Computes the solutions of a system by the Gaussian Elimination method.
    :param matrix0: A Matrix object; augmented linear system.
    :return xMatrix: A matrix object; solutions to the linear system.
    '''
    geMatrix = Matrix()
    for i in matrix0.matrix:
        row = []
        for j in i:
            row.append(j)
        geMatrix.addRow(row)
    geMatrix.columns = matrix0.columns
    geMatrix.rows = matrix0.rows

    for j in range(geMatrix.rows):
        p = 0
        greatest = 0
        for i in range(j, geMatrix.rows):   # compute pivot index
            if abs(geMatrix.matrix[i][j]) > greatest:
                greatest = geMatrix.matrix[i][j]
                p = i
        if geMatrix.matrix[p][j] == 0:
            break
        if p > j:   # interchange row p and j
            geMatrix.matrix[p], geMatrix.matrix[j] = geMatrix.matrix[j], geMatrix.matrix[p]

        for i in range(j, geMatrix.rows):   # For each i > j, subtract C[i][j]/C[j][j] times row j from row i
            eliminator = geMatrix.matrix[i][j] / geMatrix.matrix[j][j]
            if i > j:
                for x in range(geMatrix.columns):
                    geMatrix.matrix[i][x] -= eliminator * geMatrix.matrix[j][x]

    d = []
    for a in range(geMatrix.rows):  # partition C as D(n x n) and e(n x 1)
        row = []
        for b in range(geMatrix.rows):
            row.append(geMatrix.matrix[a][b])
        d.append(row)
    e = []
    for a in geMatrix.matrix:
        row = a[-1]
        e.append(row)

    x = []      # compute final solutions
    for y in range(geMatrix.rows):
        x.append(0)
    for j in range(geMatrix.rows - 1, -1, -1):
        ans = 0
        for i in range(geMatrix.rows):
            ans += d[j][i] * x[i]
        x[j] = (1/d[j][j])*(e[j] - ans)

    xMatrix = Matrix()
    for i in range(len(x)):
        row = []
        row.append(x[i])
        xMatrix.addRow(row)
        xMatrix.rows += 1
    xMatrix.columns = 1
    return xMatrix


#### Absolute Value Matrix ####
def absValMatrix(matrix0):
    ''' Creates a matrix populated with absolute values of all elements in a given matrix
    :param matrix0: A Matrix object
    :return avMatrix: A Matrix object
    '''
    avMatrix = Matrix()
    for i in matrix0.matrix:
        row = []
        for j in i:
            row.append(abs(round(j, 10)))
        avMatrix.addRow(row)
    return avMatrix


#### Magnitude of Matrix ####
def magnitude(matrix0):
    ''' Computes the magnitude of a given matrix
    :param matrix0: A Matrix object
    :return max(rowSums): A numerical data type; float/int
    '''
    rowSums = []
    for i in matrix0.matrix:
        sum = 0
        for j in i:
            sum += j
        rowSums.append(sum)
    return max(rowSums)