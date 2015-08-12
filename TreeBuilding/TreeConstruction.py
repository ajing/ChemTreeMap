"Borrow the code from GSoC2013 of Biopython"

"""Classes and methods for tree construction"""

import copy
import numpy as np
import math
import datetime
#from DataStruct import Tree
from ete2 import Tree


class Matrix(object):
    """A base class for distance matrix or scoring matrix that accepts
    a list of names and a lower triangular matrix.

    matrix = [[0],
              [1, 0],
              [2, 3, 0],
              [4, 5, 6, 0]]
    represents the symmetric matrix of
    [0,1,2,4]
    [1,0,3,5]
    [2,3,0,6]
    [4,5,6,0]

    :Parameters:
        names : list
            names of elements, used for indexing
        matrix : list
            nested list of numerical lists in lower triangular format

    Example
    -------

    >>> from Bio.Phylo.TreeConstruction import Matrix
    >>> names = ['Alpha', 'Beta', 'Gamma', 'Delta']
    >>> matrix = [[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]]
    >>> m = Matrix(names, matrix)
    >>> m
    Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [1, 0], [2, 3, 0], [4, 5, 6, 0]])

    You can use two indices to get or assign an element in the matrix.

    >>> m[1,2]
    3
    >>> m['Beta','Gamma']
    3
    >>> m['Beta','Gamma'] = 4
    >>> m['Beta','Gamma']
    4

    Further more, you can use one index to get or assign a list of elements related to that index.

    >>> m[0]
    [0, 1, 2, 4]
    >>> m['Alpha']
    [0, 1, 2, 4]
    >>> m['Alpha'] = [0, 7, 8, 9]
    >>> m[0]
    [0, 7, 8, 9]
    >>> m[0,1]
    7

    Also you can delete or insert a column&row of elemets by index.

    >>> m
    Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]])
    >>> del m['Alpha']
    >>> m
    Matrix(names=['Beta', 'Gamma', 'Delta'], matrix=[[0], [4, 0], [5, 6, 0]])
    >>> m.insert('Alpha', [0, 7, 8, 9] , 0)
    >>> m
    Matrix(names=['Alpha', 'Beta', 'Gamma', 'Delta'], matrix=[[0], [7, 0], [8, 4, 0], [9, 5, 6, 0]])

    """

    def __init__(self, names, matrix=None):
        """Initialize matrix by a list of names and a list of
        lower triangular matrix data"""
        # check names
        if isinstance(names, list) and all(isinstance(s, str) for s in names):
            if len(set(names)) == len(names):
                self.names = names
            else:
                print names
                raise ValueError("Duplicate names found")
        else:
            raise TypeError("'names' should be a list of strings")

        # check matrix
        if matrix is None:
            # create a new one with 0 if matrix is not assigned
            matrix = [[0]*i for i in range(1, len(self) + 1)]
            self.matrix = matrix
        else:
            # check if all elements are numbers
            self.matrix = [[0]*i for i in range(1, len(self) + 1)]
            for i in range(len(self)):
                for j in range(i):
                    self.matrix[i][j] = matrix[i,j]


    def __getitem__(self, item):
        """Access value(s) by the index(s) or name(s).
        For a Matrix object 'dm':
        dm[i]                   get a value list from the given 'i' to others;
        dm[i, j]                get the value between 'i' and 'j';
        dm['name']              map name to index first
        dm['name1', 'name2']    map name to index first
        """
        # Handle single indexing
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            return [self.matrix[index][i] for i in range(0, index)] + [self.matrix[i][index] for i in range(index, len(self))]
        # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            if row_index > col_index:
                return self.matrix[row_index][col_index]
            else:
                return self.matrix[col_index][row_index]
        else:
            raise TypeError("Invalid index type.")

    def __setitem__(self, item, value):
        """Set value by the index(s) or name(s).
        Similar to __getitem__
        dm[1] = [1, 0, 3, 4]    set values from '1' to others;
        dm[i, j] = 2            set the value from 'i' to 'j'
        """

        # Handle single indexing
        if isinstance(item, (int, str)):
            index = None
            if isinstance(item, int):
                index = item
            elif isinstance(item, str):
                if item in self.names:
                    index = self.names.index(item)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if index > len(self) - 1:
                raise IndexError("Index out of range.")
            # check and assign value
            if isinstance(value, list) and all(isinstance(n, (int, long, float, complex)) for n in value):
                if len(value) == len(self):
                    for i in range(0, index):
                        self.matrix[index][i] = value[i]
                    for i in range(index, len(self)):
                        self.matrix[i][index] = value[i]
                else:
                    raise ValueError("Value not the same size.")
            else:
                raise TypeError("Invalid value type.")
        # Handle double indexing
        elif len(item) == 2:
            row_index = None
            col_index = None
            if all(isinstance(i, int) for i in item):
                row_index, col_index = item
            elif all(isinstance(i, str) for i in item):
                row_name, col_name = item
                if row_name in self.names and col_name in self.names:
                    row_index = self.names.index(row_name)
                    col_index = self.names.index(col_name)
                else:
                    raise ValueError("Item not found.")
            else:
                raise TypeError("Invalid index type.")
            # check index
            if row_index > len(self) - 1 or col_index > len(self) - 1:
                raise IndexError("Index out of range.")
            # check and assign value
            if isinstance(value, (int, long, float, complex)):
                if row_index > col_index:
                    self.matrix[row_index][col_index] = value
                else:
                    self.matrix[col_index][row_index] = value
            else:
                raise TypeError("Invalid value type.")
        else:
            raise TypeError("Invalid index type.")

    def __delitem__(self, item):
        """Delete related distances by the index or name"""
        index = None
        if isinstance(item, int):
            index = item
        elif isinstance(item, str):
            index = self.names.index(item)
        else:
            raise TypeError("Invalid index type.")
        # remove distances related to index
        for i in range(index + 1, len(self)):
            del self.matrix[i][index]
        del self.matrix[index]
        # remove name
        del self.names[index]

    def insert(self, name, value, index=None):
        """Insert distances given the name and value.

        :Parameters:
        name : str
            name of a row/col to be inserted
        value : list
            a row/col of values to be inserted
        """
        if isinstance(name, str):
            # insert at the given index or at the end
            if index is None:
                index = len(self)
            if not isinstance(index, int):
                raise TypeError("Invalid index type.")
            # insert name
            self.names.insert(index, name)
            # insert elements of 0, to be assigned
            self.matrix.insert(index, [0] * index)
            for i in range(index, len(self)):
                self.matrix[i].insert(index, 0)
            # assign value
            self[index] = value
        else:
            raise TypeError("Invalid name type.")

    def __len__(self):
        """Matrix length"""
        return len(self.names)

    def __repr__(self):
        return self.__class__.__name__ \
        + "(names=%s, matrix=%s)" \
        % tuple(map(repr, (self.names, self.matrix)))

    def __str__(self):
        """Get a lower triangular matrix string"""
        matrix_string = '\n'.join([self.names[i] + "\t" +
            "\t".join([str(n) for n in self.matrix[i]]) for i in range(0, len(self))])
        matrix_string = matrix_string + "\n\t" + "\t".join(self.names)
        return matrix_string


class DistanceMatrix(Matrix):
    """Distance matrix class that can be used for distance based tree
     algorithms.
    All diagonal elements will be zero no matter what the users provide.
    """

    def __init__(self, names, matrix=None):
        Matrix.__init__(self, names, matrix)
        self._set_zero_diagonal()

    def __setitem__(self, item, value):
        Matrix.__setitem__(self, item, value)
        self._set_zero_diagonal()

    def _set_zero_diagonal(self):
        """set all diagonal elements to zero"""
        for i in range(0, len(self)):
            self.matrix[i][i] = 0

def AddLineForNode(clade, moldict, filecontent):
    node_size  = None
    node_color = None
    alpha      = 0.3
    #print moldict[clade.name]
    try:
        node_size = str(math.log(moldict[clade.name][0] + 1, 100) * alpha)
        if moldict[clade.name][1] in ["allosteric", "competitive"]:
            node_color = moldict[clade.name][1] == "allosteric" and "red" or "blue"
        else:
            #node_color = str(moldict[clade.name][1])
            node_color = str(moldict[clade.name][2])
        node_line = clade.name + "[label=\"\", width=" + node_size + " color=" + node_color + " ];"
    except:
        node_line = clade.name + "[label=\"\", width=0 ];"
    filecontent.append(node_line)

def AddRelation(cladeChild, cladeParent, filecontent):
    node1_relation = cladeChild.name + " -- " + cladeParent.name + " [len=" + "{:f}".format(cladeChild.dist).rstrip("0") + "]"
    #node1_relation = cladeChild.name + " -- " + cladeParent.name
    filecontent.append(node1_relation + ";")

def AddTwoChild(cladeChild1, cladeChild2, cladeParent, filecontent):
    AddRelation(cladeChild1, cladeParent, filecontent)
    AddRelation(cladeChild2, cladeParent, filecontent)

def nj(distance_matrix, moldict, outfilename = False):
    """Construct and return an Neighbor Joining tree.

    :Parameters:
        distance_matrix : DistanceMatrix
            The distance matrix for tree construction.
    """

    # file content
    filecontent = []
    aline = "graph G{\nnode [shape=circle, style=filled];"
    filecontent.append(aline)

    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("Must provide a DistanceMatrix object.")

    # make a copy of the distance matrix to be used
    dm = copy.deepcopy(distance_matrix)
    # init terminal clades
    clades = [Tree(name = name) for name in dm.names]
    # init node distance
    node_dist = [0] * len(dm)
    # init minimum index
    min_i = 0
    min_j = 0
    times = 0
    while len(dm) > 2:
        times += 1
        print "dimension of Distance Matrix:", len(dm)
        #if len(dm) < 8:
        #    print dm
        # calculate nodeDist
        for i in range(0, len(dm)):
            node_dist[i] = 0
            for j in range(0, len(dm)):
                node_dist[i] += dm[i, j]
            node_dist[i] = node_dist[i] / (len(dm) - 2)

        # find minimum distance pair
        min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
        min_i = 0
        min_j = 1
        for i in range(1, len(dm)):
            for j in range(0, i):
                temp = dm[i, j] - node_dist[i] - node_dist[j]
                if min_dist > temp:
                    min_dist = temp
                    min_i = i
                    min_j = j
        # create clade
        clade1 = clades[min_i]
        clade2 = clades[min_j]
        inner_clade = Tree(name = "_".join([clade1.name, clade2.name]))
        inner_clade.add_child(clade1)
        inner_clade.add_child(clade2)
        #assign branch length
        clade1.dist = (dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]) / 2
        clade2.dist = dm[min_i, min_j] - clade1.dist

        # add lines to dot file
        if not "_" in clade1.name:
            AddLineForNode(clade1, moldict, filecontent)
        if not "_" in clade2.name:
            AddLineForNode(clade2, moldict, filecontent)
        # relationship between nodes
        if len(dm) > 3:
            AddLineForNode(inner_clade, moldict, filecontent)
            AddTwoChild(clade1, clade2, inner_clade, filecontent)

        # update node list
        clades[min_j] = inner_clade
        del clades[min_i]

        # rebuild distance matrix,
        # set the distances of new node at the index of min_j
        for k in range(0, len(dm)):
            if k != min_i and k != min_j:
                dm[min_j, k] = (dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]) * 1.0 / 2

        dm.names[min_j] = "_".join([clade1.name, clade2.name])

        del dm[min_i]

    # set the last clade as one of the child of the inner_clade
    root = None
    if clades[0] == inner_clade:
        clades[0].dist = 0
        clades[1].dist = dm[1, 0]
        clades[0].add_child(clades[1])
        clades[0].name = clades[0].name + "_" + clades[1].name
        AddLineForNode(clades[0], moldict, filecontent)
        AddLineForNode(clades[1], moldict, filecontent)
        AddTwoChild(clade1, clade2, clades[0], filecontent)
        AddRelation(clades[1], clades[0], filecontent)
        root = clades[0]
    else:
        clades[0].dist = dm[1, 0]
        clades[1].dist = 0
        clades[1].add_child(clades[0])
        clades[1].name = clades[1].name + "_" + clades[0].name
        AddLineForNode(clades[0], moldict, filecontent)
        AddLineForNode(clades[1], moldict, filecontent)
        AddTwoChild(clade1, clade2, clades[1], filecontent)
        AddRelation(clades[0], clades[1], filecontent)
        root = clades[1]

    # create dot langauge to file:
    fmt='./Data/%Y-%m-%d-%Hh-%Mm_dot.gv'
    newfilename = datetime.datetime.now().strftime(fmt)
    newfileobj  = open(newfilename, "w")
    # close file obj for dot language
    newfileobj.write("\n".join(filecontent) + "}")
    newfileobj.close()

    if outfilename:
        return newfilename
    else:
        return root

if __name__ == "__main__":
    # some test cases
    #matrixfile = "./Data/similarityMatrix_small.npy"
    #dmatrix    = np.load(matrixfile)
    dmatrix    = np.matrix("1 0 0 0; 14 1 0 0 ; 12 18 1 0 ; 17 21 27 1")
    dmatrix    = np.matrix("1 0 0 0; 17 1 0 0 ; 21 12 1 0 ; 27 18 14 1")
    #dlist      = ["1","2","3","4","5","6","7","8","9"]
    dlist      = ["A","B","C","D"]
    moldict    = {"A":[4,"allosteric"],"B":[3,"competitive"],"C":[6,"allosteric"],"D":[1,"competitive"]}
    DMatrix    = DistanceMatrix(dlist, dmatrix)
    print nj(DMatrix, moldict).write(format = 7)
    print nj(DMatrix, moldict)
