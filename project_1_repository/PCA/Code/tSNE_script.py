import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

def tSNE_impl(file_name):
    #read the last column of the file as disease list
    file_handle = open(file_name + ".txt", "r")
    lines = file_handle.readlines()
    num_cols = len(lines[0].strip().split("\t"))

    disease_list = []

    for line in lines:
        disease_list.append(line.strip().split("\t")[-1])

    #encode the disease list assigning a unique number to each disease in order to color the final plot
    d = dict()
    counter = 0
    for disease in disease_list:
        if disease in d:
            continue
        d[disease] = counter
        counter += 1
    disease_list_encoded = [d[disease] for disease in disease_list]


    #get list of unique diseases
    unique_disease_set = set(disease_list)
    unique_disease_list = list(unique_disease_set)

    #get list of encodings
    unique_disease_encoded_set = set(disease_list_encoded)
    unique_disease_encoded_list = list(unique_disease_encoded_set)

    #read the input features into numpy array
    input_matrix = np.loadtxt(file_name + ".txt", delimiter="\t", usecols = range(num_cols - 1))

    tsne = TSNE(n_components=2, n_iter=600)
    tsne_results = tsne.fit_transform(input_matrix)

    # plot diseases using the reduced matrix as the coordinates
    colors = [plt.cm.jet(float(i) / max(unique_disease_encoded_list)) for i in unique_disease_encoded_list]
    for i, u in enumerate(unique_disease_list):
        xi = [p for (j,p) in enumerate(tsne_results[:,0]) if disease_list[j] == u]
        yi = [p for (j,p) in enumerate(tsne_results[:,1]) if disease_list[j] == u]
        plt.scatter(xi, yi, c=colors[i], label=str(u))

    plt.title(file_name + " tSNE scatter plot")
    plt.legend()
    plt.show()

#tSNE_impl("pca_a")
#tSNE_impl("pca_b")
#tSNE_impl("pca_c")
tSNE_impl("pca_demo")