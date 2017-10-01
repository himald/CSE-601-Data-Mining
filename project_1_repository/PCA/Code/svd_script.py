import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg;

def svd_impl(file_name):
    #read the last column of the file as disease list
    file_handle = open(file_name + ".txt", "r")
    lines = file_handle.readlines()
    num_cols = len(lines[0].strip().split("\t"))

    disease_list = []

    for line in lines:
        disease_list.append(line.strip().split("\t")[-1])

    # encode the disease list assigning a unique number to each disease in order to color the final plot
    d = dict()
    counter = 0
    for disease in disease_list:
        if disease in d:
            continue
        d[disease] = counter
        counter += 1
    disease_list_encoded = [d[disease] for disease in disease_list]


    # get list of unique diseases
    unique_disease_set = set(disease_list)
    unique_disease_list = list(unique_disease_set)

    # get list of encodings
    unique_disease_encoded_set = set(disease_list_encoded)
    unique_disease_encoded_list = list(unique_disease_encoded_set)

    # read the input features into numpy array
    input_matrix = np.loadtxt(file_name + ".txt", delimiter="\t", usecols = range(num_cols - 1))

    # obtain U, s, v matrices by svd of input matrix
    U, s, V = linalg.svd(input_matrix)

    # obtain reduced matrix by selecting first n columns in matrix U
    dimensions = 2
    reduced_matrix = U[:,:dimensions]

    # plot diseases using the reduced matrix as the coordinates
    colors = [plt.cm.jet(float(i) / max(unique_disease_encoded_list)) for i in unique_disease_encoded_list]
    for i, u in enumerate(unique_disease_list):
        xi = [p for (j,p) in enumerate(reduced_matrix[:,0]) if disease_list[j] == u]
        yi = [p for (j,p) in enumerate(reduced_matrix[:,1]) if disease_list[j] == u]
        plt.scatter(xi, yi, c=colors[i], label=str(u))

    plt.title(file_name + " SVD scatter plot")
    plt.legend()
    plt.show()

#svd_impl("pca_a")
#svd_impl("pca_b")
#svd_impl("pca_c")
svd_impl("pca_demo")