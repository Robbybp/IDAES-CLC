import numpy as np 
import scipy.linalg as la
import matplotlib.pyplot as plt
import csv
import os

def make_dir(filename):
    if not os.path.isdir('profiles'):
        os.mkdir('profiles')
    with open(filename, newline='') as f:
        reader = csv.reader(f, delimiter=',')
        header = next(reader)
        n_idx = 0
        for h in header:
            if 'idx' in h:
                n_idx = n_idx + 1
            if h.isnumeric():
                break
        print(n_idx)
        for row in reader:
            name = row[0]
            idx_list = []
            for idx in row[1:n_idx+1]:
                if idx != '':
                    idx_list.append(idx)
            if len(idx_list) == 1:
                idx_label = idx_list[0]
            elif len(idx_list) > 1:
                idx_label = idx_list[0]
                for idx in idx_list[1:len(idx_list)]:
                    idx_label = idx_label + '_' + idx

            var_rel_path = 'profiles/' + name
            if not os.path.isdir(var_rel_path):
                os.mkdir(var_rel_path)
            idx_rel_path = var_rel_path + '/' + idx_label
            if not os.path.isdir(idx_rel_path):
                os.mkdir(idx_rel_path)
        
def get_fig_filename(filename,row,n_index):
    filename = filename[0:len(filename)-4]
    # ^ extract the non-'.csv' part of the filename
    name = row[0]
    idx_list = []
    for idx in row[1:n_index+1]:
        if idx != '':
            idx_list.append(idx)
    if len(idx_list) == 1:
        idx_label = idx_list[0]
    elif len(idx_list) > 1:
        idx_label = idx_list[0]
        for idx in idx_list[1:len(idx_list)]:
            idx_label = idx_label + '_' + idx

    fig_rel_path = 'profiles/' + name + '/' + idx_label + '/' + filename + '.png'
    return fig_rel_path

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_data(row, n_index):
    data = []
    for entry in row[n_index+1:]:
        data.append(float(entry))
    return data

def plot_data(data, time, filename=None):
    fig,ax = plt.subplots()
    ax.plot(time,data)
    if filename != None:
        fig.savefig(filename)
        print('Saved fig at ' + filename)
    return fig
        
def construct_plots(filename):
    # input is the csv file name
    make_dir(filename)
    with open(filename, newline='') as f:
        reader = csv.reader(f, delimiter=',')
        header = next(reader)
        n_idx = 0
        for h in header:
            if 'idx' in h:
                n_idx = n_idx + 1
            if h.isnumeric():
                break
        time = get_data(header, n_idx)
        for row in reader:
            data = get_data(row, n_idx)
            plt_filename = get_fig_filename(filename, row, n_idx)
            plot_data(data, time, plt_filename)
            plt.cla()
            plt.clf()


def main():
    print('executing main...')
    filename = 'gas_T.csv'
    construct_plots(filename)

if __name__ == '__main__':
    main()
