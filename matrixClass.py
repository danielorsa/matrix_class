# This module employs a class for creating matrices and performing some matrix operations

class Matrix:
    def __init__(self):
        self.totalElements = 0
        self.rows = 0
        self.columns = 0
        self.matrix = []
        self.classNum = None

    def setClassNum(self, cn):
        self.classNum = cn

    def intake(self, file):
        ''' Populates a matrix with data from a text file
        :param file: A text file; string
        :return none
        '''
        infile = open(file, "r")
        for line in infile:
            row = []
            elements = line.split()
            for i in elements:
                if "." in i:
                    row.append(float(i))
                else:
                    row.append(int(i))
                self.totalElements += 1
            self.matrix.append(row)
            self.rows += 1
        self.columns = int(self.totalElements / self.rows)

    def intakeFromList(self, dataList):
        self.matrix = dataList
        self.rows = len(dataList)
        self.columns = len(dataList[0])

    def identity(self, n):
        ''' Creates an identity matrix; Only call after creating a new empty matrix
        :param n: dimension of identity matrix; integer
        :return none
        '''
        self.rows = n
        self.columns = n
        self.totalElements = n * n
        pos = 0
        for i in range(n):
            row = []
            for j in range(n):
                row.append(0)
            self.matrix.append(row)
            self.matrix[pos][pos] = 1
            pos += 1

    def getMatrix(self):
        ''' Prints contents of matrix
        :param none
        :return none
        '''
        for line in self.matrix:
            print(line)

    def addRow(self, row):
        ''' Adds new row to end of matrix
        :param none
        :return none
        '''
        self.matrix.append(row)

    def inverse(self):
        ''' Computes and returns inverse of a matrix by Gauss-Jordan method
        :param none
        :return invMatrix: a Matrix object
        '''
        idMatrix = Matrix()     # create identity matrix
        idMatrix.identity(self.columns)
        invMatrix = Matrix()
        invMatrix.rows = self.rows
        invMatrix.columns = self.columns * 2

        for line in self.matrix:
            invMatrix.addRow(line)
        for i in range(self.rows):  # augment coefficient matrix with identity matrix
            for j in range(idMatrix.columns):
                invMatrix.matrix[i].append(idMatrix.matrix[i][j])

        for j in range(invMatrix.rows):
            p = 0
            greatest = 0
            for i in range(j, invMatrix.rows):  # compute pivot index
                if abs(invMatrix.matrix[i][j]) > greatest:
                    greatest = invMatrix.matrix[i][j]
                    p = i
            if invMatrix.matrix[p][j] == 0:
                break
            if p > j:   # interchange row p and j
                invMatrix.matrix[p], invMatrix.matrix[j] = invMatrix.matrix[j], invMatrix.matrix[p]
            for i in range(invMatrix.columns):  # divide row j by C[j][j]
                newElement = invMatrix.matrix[j][i] / greatest
                invMatrix.matrix[j][i] = newElement
            for i in range(invMatrix.rows):     # for each i =/= j, subtract C[i][j] times row j from row i
                eliminator = invMatrix.matrix[i][j]
                if i != j:
                    for x in range(invMatrix.columns):
                        invMatrix.matrix[i][x] -= eliminator * invMatrix.matrix[j][x]

        for i in range(invMatrix.rows):     # partition inverse matrix
            newRow = []
            for j in range(invMatrix.rows):
                newRow.append(invMatrix.matrix[i][j+invMatrix.rows])
            invMatrix.matrix[i] = newRow

        invMatrix.columns = self.columns

        return invMatrix


    def determinant(self):
        ''' Computes and returns determinant of a matrix by Gaussian elimination method
        :param none
        :return determinant: A numerical data type; int/float
        '''
        detMatrix = Matrix()
        detMatrix.rows = self.rows
        detMatrix.columns = self.columns
        for line in self.matrix:
            detMatrix.addRow(line)

        r = 0
        for j in range(detMatrix.rows):
            p = 0
            greatest = 0
            for i in range(j, detMatrix.rows):      # compute pivot index
                if abs(detMatrix.matrix[i][j]) > greatest:
                    greatest = detMatrix.matrix[i][j]
                    p = i
            if detMatrix.matrix[p][j] == 0:
                determ = 0
                break
            if p > j:       # interchange row p and j, increment r
                detMatrix.matrix[p], detMatrix.matrix[j] = detMatrix.matrix[j], detMatrix.matrix[p]
                r += 1

            for i in range(j, detMatrix.rows):      # For each i > j, subtract A[i][j]/A[j][j] times row j from row i

                eliminator = detMatrix.matrix[i][j] / detMatrix.matrix[j][j]
                if i > j:
                    for x in range(detMatrix.columns):
                        detMatrix.matrix[i][x] -= eliminator * detMatrix.matrix[j][x]

        sign = (-1) ** r
        determ = 1
        for j in range(detMatrix.rows):     # multiply product of diagonal values by (-1)^(number of row interchanges)
            determ *= detMatrix.matrix[j][j]
        determ *= sign
        return determ

    def __str__(self):
        return "{0}x{1} Matrix".format(self.rows, self.columns)