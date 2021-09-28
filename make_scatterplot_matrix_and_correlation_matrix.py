from linear_algebra import shape, get_row, get_column, make_matrix, \
    vector_mean, vector_sum, dot, magnitude, vector_subtract, scalar_multiply
from stats import correlation, standard_deviation, mean
import matplotlib.pyplot as plt
import math, random, csv, glob

def get_subject_data(path):

    data = []

    # glob.glob returns every filename that matches the wildcarded path
    for fn in glob.glob(path): 

        with open(fn,'r') as file:
            for line in file:
                line = line.rstrip('\n')
                data.append([int(x) for x in line.split('\t')])
    return data


def correlation_matrix(data):
    """returns the num_columns x num_columns matrix whose (i, j)th entry
    is the correlation between columns i and j of data"""

    _, num_columns = shape(data)
 

    def matrix_entry(i, j):
        return correlation(get_column(data, i), get_column(data, j))

    return make_matrix(num_columns, num_columns, matrix_entry)

def make_scatterplot_matrix(path):
    data = get_subject_data(path)      
    _, num_columns = shape(data)
    fig, ax = plt.subplots(num_columns, num_columns)

    for i in range(num_columns):
        for j in range(num_columns):

            # scatter column_j on the x-axis vs column_i on the y-axis
            if i != j: ax[i][j].scatter(get_column(data, j), get_column(data, i))

            # unless i == j, in which case show the series name
            else: ax[i][j].annotate("extraction " + str(i), (0.5, 0.5),
                                    xycoords='axes fraction',
                                    ha="center", va="center")

            # then hide axis labels except left and bottom charts
            if i < num_columns - 1: ax[i][j].xaxis.set_visible(False)
            if j > 0: ax[i][j].yaxis.set_visible(False)

    # fix the bottom right and top left axis labels, which are wrong because
    # their charts only have text in them
    ax[-1][-1].set_xlim(ax[0][-1].get_xlim())
    plt.savefig('CD4_scatterplot_matrix.pdf')
    
    
    with open('CD4_correlation_matrix.txt', 'w') as f:
        print(correlation_matrix(data), file=f)

if __name__ == "__main__":

    make_scatterplot_matrix(r"/Volumes/Seagate/Mellors/DeepDive_processing_reanalysis_2020/Modeling/correlation_matrix/CN/CD4_clone_size_by_extractions_for_correlation_matrix.txt")
