import numpy as np

# Query access matrix - use of attributes in application queries
# Q1: Select BUDGET from PROJ where PNO = VALUE;
# Q2: Select PNAME, BUDGET from PROJ;
# Q3: Select PNAME from PROJ where LOC = VALUE
# Q4: Select SUM(BUDGET) from PROJ where LOC = VALUE
# A1 = PNO, A2 = PNAME. A3 = BUDGET, A4 = LOC


# query_access_matrix = [[0, 1, 1, 0, 1],
#                        [1, 1, 1, 0, 1],
#                        [1, 0, 0, 1, 1],
#                        [0, 0, 1, 0, 0],
#                        [1, 1, 1, 0, 0]]
query_access_matrix = [[1, 0, 1, 0],
                       [0, 1, 1, 0],
                       [0, 1, 0, 1],
                       [0, 0, 1, 1]]
query_access_matrix = np.array(query_access_matrix)
print("Query Access Matrix = ")
print(query_access_matrix)

# Frequency Access Matrix - no of times the queries accesses the sites in a day
# Assuming random values
# Frequency_access_matrix = [[10, 20, 0],
#                            [5, 0, 10],
#                            [0, 35, 5],
#                            [0, 10, 0],
#                            [0, 15, 0]]

Frequency_access_matrix = [[15, 20, 10],  # Q1
                           [5, 0, 0],  # Q2
                           [25, 25, 25],  # Q3
                           [3, 0, 0]]  # q4
Frequency_access_matrix = np.array(Frequency_access_matrix)
print("Frequency Access Matrix = ")
print(Frequency_access_matrix)

# Taking the sum of attr access by each query from the frequency access matrix
sum_attr_access = np.sum(Frequency_access_matrix, axis=1)
print("Sum of attr access by each query = ")
print(sum_attr_access)

no_of_queries = query_access_matrix.shape[0]  # no of queries
no_of_attr = query_access_matrix.shape[1]  # no of attributes
attr_affinity_matrix = np.zeros((no_of_attr + 1, no_of_attr))  # Attribute affinity matrix

# Setting of attribute number in the first row of the attribute affinity matrix
for i in range(no_of_attr):
    attr_affinity_matrix[0, i] = i + 1
print("attribute affinity matrix = ")
print(attr_affinity_matrix)


# Finding out the attribute usage with frequency
# called the Attribute Affinity matrix (AA Matrix)
def affinity_calc():
    global attr_affinity_matrix
    for col_attr in range(no_of_attr):
        for row_attr in range(1, no_of_attr + 1):
            affinity_value = 0
            for q in range(no_of_queries):
                if query_access_matrix[q][col_attr] == 1 & query_access_matrix[q][row_attr - 1] == 1:
                    affinity_value += sum_attr_access[q]
            attr_affinity_matrix[row_attr][col_attr] = affinity_value
    return attr_affinity_matrix


attr_affinity_matrix = affinity_calc()
print("attribute affinity matrix = ")
print(attr_affinity_matrix)


# Function to calculate the bond between two columns
def bond(left, right):
    bond_value = 0
    # Boundary conditions
    if left == -1 or left == no_of_attr or right == -1 or right == no_of_attr:
        return bond_value
    else:
        bond_value = np.sum(np.multiply(attr_affinity_matrix[1:, left], attr_affinity_matrix[1:, right]))
        return bond_value


# Function to calculate the contribution of a certain configuration of columns.
# Example: to calculate contribution of placement of (A2, A4, A3), this function is called with
# left = 2, middle = 4, right = 3, and a reference to the Affinity Matrix object.
def contribution(left, middle, right):
    a = bond(left, middle)
    b = bond(middle, right)
    c = bond(left, right)
    if right == middle + 1:
        cont = 2 * a
    else:
        cont = 2 * (a + b - c)
    print("cont( A", left + 1, ", A", middle + 1, ", A", right + 1, ") = 2 * (", a, " + ", b,
          " - ", c, ") = ", cont)
    return cont


# Bond Energy Algorithm
def bea_algo():
    clustered_affinity_matrix = np.zeros((no_of_attr + 1, no_of_attr))
    # Copy the first and second columns from the Attribute Affinity Matrix
    clustered_affinity_matrix[:, 0] = attr_affinity_matrix[:, 0]
    clustered_affinity_matrix[:, 1] = attr_affinity_matrix[:, 1]
    print("After shifting rows 0 and 1=")
    print(clustered_affinity_matrix)

    # Best Placement of attributes starting from index 2 that is the third attribute
    for index in range(2, no_of_attr):
        contribution_array = []
        print("Best location for attribute = ", attr_affinity_matrix[:, index])
        # Calculating the contribution value
        for i in range(index):
            contribution_value = contribution(i - 1, index, i)
            contribution_array.append(contribution_value)

        contribution_array.append(contribution(index - 1, index, index + 1))
        # loc <- placement given by max contribution value
        loc = contribution_array.index(max(contribution_array))
        print("Location of max cont = ", loc + 1)

        # Shifting attribute to the location of max contribution in CA
        for k in range(index, loc, -1):
            clustered_affinity_matrix[:, k] = clustered_affinity_matrix[:, k - 1]
        clustered_affinity_matrix[:, loc] = attr_affinity_matrix[:, index]
        print("CA after swapping attribute", index + 1)
        print(clustered_affinity_matrix)

        # Shifting attribute to the location of max contribution in AA
        temp = attr_affinity_matrix[:, index].copy()
        for m in range(index, loc, -1):
            attr_affinity_matrix[:, m] = attr_affinity_matrix[:, m - 1]
        attr_affinity_matrix[:, loc] = temp
        print("AA after swapping attribute", index + 1)
        print(attr_affinity_matrix)

    # Interchanging of rows in CA after the BEA algorithm
    CA_ordered_row = np.zeros((no_of_attr, no_of_attr))
    n = 0
    for m in range(no_of_attr):
        order = clustered_affinity_matrix[0, :]
        CA_ordered_row[n, :] = clustered_affinity_matrix[int(order[m]), :]
        n += 1
    clustered_affinity_matrix[1:][:] = CA_ordered_row
    print("CA after interchanging rows = ")
    print(clustered_affinity_matrix)

    return clustered_affinity_matrix


CA = bea_algo()
print("CA =")
print(CA)
