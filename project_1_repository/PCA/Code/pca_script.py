import numpy as np
import matplotlib.pyplot as plt


def pca(file_name):
    # read the last column of the file as disease list
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

    # find adjusted_matrix
    adjusted_matrix = input_matrix - input_matrix.mean(axis=0)

    # covariance
    # cov_matrix = np.dot(adjusted_matrix.T, adjusted_matrix) / adjusted_matrix.shape[0]
    cov_matrix = np.cov(adjusted_matrix.T)

    # obtain eigen values and eigen vectors of covariance matrix
    eig_val, eig_vec = np.linalg.eig(cov_matrix)

    # select top n eigen values as the principal components
    top_2_eig_val_indexes = eig_val.argsort()[-2:][::-1]
    top_2_eig_vec = eig_vec[:, top_2_eig_val_indexes]
    principle_components_matrix = np.empty([input_matrix.shape[0], top_2_eig_vec.shape[1]])
    i = 0
    for top_eig_vec in top_2_eig_vec.T:
        principle_components_matrix[:,i] = np.dot(adjusted_matrix, top_eig_vec.T)
        i += 1

    # plot diseases using the principle components as the coordinates
    colors = [plt.cm.jet(float(i) / max(unique_disease_encoded_list)) for i in unique_disease_encoded_list]
    for i, u in enumerate(unique_disease_list):
        xi = [p for (j,p) in enumerate(principle_components_matrix[:,0]) if disease_list[j] == u]
        yi = [p for (j,p) in enumerate(principle_components_matrix[:,1]) if disease_list[j] == u]
        plt.scatter(xi, yi, c=colors[i], label=str(u))

    plt.title(file_name + " scatter plot")
    plt.legend()
    plt.show()

#pca("pca_a")
#pca("pca_b")
#pca("pca_c")
pca("pca_demo")

