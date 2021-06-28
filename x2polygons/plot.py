""" This module contains all plot related functionality.
"""

import matplotlib.pyplot as plt
import subprocess, os

def plot_turn_function(turn):  
    '''
    Plots the cumulative length turn function. 
    
    Args:
        - **turn** (*dict*): The turn dictionary of a polygon.
    
    Returns:
        - Plotted turn function 
    '''   
    # Plot the turn function
    cum_sum_lengths = [] + turn["lengths"] 
    cum_sum_angles = [] + turn["angles"]
        
    for i in range(1,len(turn["lengths"])):
        cum_sum_angles[i] = cum_sum_angles[i] + cum_sum_angles[i-1]
        cum_sum_lengths[i] = cum_sum_lengths[i] + cum_sum_lengths[i-1]
    
    #multiply with two to obtain the piecewise nature to plot
    piece_wise_angles = [0]*(len(cum_sum_angles)*2)
    piece_wise_lengths = [0]*(len(cum_sum_angles)*2)
    
    for i in range(len(cum_sum_angles)):
        # y
        piece_wise_angles[i*2] = cum_sum_angles[i]
        piece_wise_angles[i*2 + 1] = cum_sum_angles[i]
    
    for i in range(len(cum_sum_lengths)-1):
        # x
        piece_wise_lengths[i*2 + 1] = cum_sum_lengths[i]  
        piece_wise_lengths[i*2 + 2] = cum_sum_lengths[i]
    piece_wise_lengths[-1] = cum_sum_lengths[-1]
    
    
    plt.plot(piece_wise_lengths, piece_wise_angles)
    plt.xlabel("Length of edges")
    plt.ylabel("Turn angles")
    plt.show()

def plot_polygon(poly_a):
    '''
    Plots an input polygon
    
    Args:
        - **poly_a** (*polygon*): The input polygon to be plotted.
    
    Returns:
        - The plotted polygon
    '''
    fig, ax = plt.subplots()
    
    # Plot the edges
    for i in range(len(poly_a.exterior.coords)-1):
        ax.plot([poly_a.exterior.coords[i][0], poly_a.exterior.coords[i+1][0]], 
                 [poly_a.exterior.coords[i][1], poly_a.exterior.coords[i+1][1]], 'b')
    
    # Plot the nodes
    for i in range(len(poly_a.exterior.coords)-1):
        ax.plot(poly_a.exterior.coords[i][0], 
                 poly_a.exterior.coords[i][1], 'bo',
                 markersize = 12,
                 fillstyle = 'none')
    
    # Remove the axes
    fig.patch.set_visible(False)
    ax.axis('off')

def plot_x2polygons(poly_a, poly_b, **kwargs):
    '''
    Plots two matching (homologous) polygons - one from OSM and the other from the reference dataset.
    
    Args:
        - **poly_a** (*polygon*): First polygon
        - **poly_b** (*polygon*): Second polygon
        - *kwargs*: 
            - **save_plot_as** (*str*): Output filename as svg
    
    Returns:
        - Plot of the polygons
        
    Examples:
        
        >>> plot_x2polygons(poly_a, poly_b, save_plot_as = "my_output.svg")
    '''
    fig, ax = plt.subplots()
    
    # First polygon
    # Plot the edges 
    for i in range(len(poly_a.exterior.coords)-1):
        ax.plot([poly_a.exterior.coords[i][0], poly_a.exterior.coords[i+1][0]], 
                 [poly_a.exterior.coords[i][1], poly_a.exterior.coords[i+1][1]], 'r')
        
    # Plot the nodes
    for i in range(len(poly_a.exterior.coords)-1):
        ax.plot(poly_a.exterior.coords[i][0], 
                 poly_a.exterior.coords[i][1], 'rs',
                 markersize = 12,
                 fillstyle = 'full',
                 label = "A")
    
    # Second polygon
    # Plot the edges 
    for i in range(len(poly_b.exterior.coords)-1):
        ax.plot([poly_b.exterior.coords[i][0], poly_b.exterior.coords[i+1][0]], 
                 [poly_b.exterior.coords[i][1], poly_b.exterior.coords[i+1][1]], 'b',
                 linestyle='dashed')
        
    # Plot the nodes
    for i in range(len(poly_b.exterior.coords)-1):
        ax.plot(poly_b.exterior.coords[i][0], 
                 poly_b.exterior.coords[i][1], 'bo',
                 markersize = 8,
                 fillstyle = 'full',
                 label = 'B')
    
    # Remove the axes
    fig.patch.set_visible(False)
    ax.axis('off')
    
    # Make sure the length of one unit in X & Y axis is the same
    ax.set_aspect('equal')
    
    # Setting the axis limits
    # ax.set_xlim(-2, 12)
    # ax.set_ylim(-2, 12)
    
    
    # Do not repeat the legend labels (i.e. for each node).
    # REF: https://stackoverflow.com/questions/13588920/stop-matplotlib-repeating-labels-in-legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc = 'center')
    
    
    # Node labels 
    save_nodes = []
    if ('with_node_labels' in kwargs):
        # Position the labels w.r.t. to the location of the nodes - how much to shift?
        label_drift = kwargs['with_node_labels']
        
        # For polygon A
        for i in range(len(poly_a.exterior.coords)-1):
            x_val = str(int(poly_a.exterior.coords[i][0]))
            y_val = str(int(poly_a.exterior.coords[i][1]))
            
            node_label = 'a' + str(i) + "(" + x_val + "," + y_val + ")"  
            ax.text(poly_a.exterior.coords[i][0]-label_drift[0], poly_a.exterior.coords[i][1]-label_drift[0], 
                    node_label, 
                    style='italic',
                    color='red')
            save_nodes.append((poly_a.exterior.coords[i][0], poly_a.exterior.coords[i][1]))
        
        # For Polygon B
        for i in range(len(poly_b.exterior.coords)-1):
            # If a node B coincides with a node from A, skip its coordinates
            x_val = str(int(poly_b.exterior.coords[i][0]))
            y_val = str(int(poly_b.exterior.coords[i][1]))
            
            node_tmp = (poly_b.exterior.coords[i][0], poly_b.exterior.coords[i][1])
            if (node_tmp in save_nodes):
                node_label = 'b' + str(i)
            else:
                node_label = 'b' + str(i) + "(" + x_val + "," + y_val + ")"  
            
            ax.text(poly_b.exterior.coords[i][0]+label_drift[1], poly_b.exterior.coords[i][1]+label_drift[1], 
                    node_label, 
                    style='italic',
                    color='blue')
        
            #ax.text -> bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10}

    # Save the plot as a vector graphic
    # REF: https://stackoverflow.com/questions/9266150/matplotlib-generating-vector-plot
    if('file_path' in kwargs):
        inkscape_path = kwargs.get('inkscape', "C://Program Files//Inkscape//bin//inkscape.exe")
        filepath = kwargs.get('file_path', None)
        
        if filepath is not None:
            path, filename = os.path.split(filepath)
            filename, extension = os.path.splitext(filename)
    
            svg_filepath = os.path.join(path, filename+'.svg')
            emf_filepath = os.path.join(path, filename+'.emf')
    
            fig.savefig(svg_filepath, format='svg')
    
            subprocess.call([inkscape_path, svg_filepath, '--export-filename', emf_filepath])
            os.remove(svg_filepath)
            
        

    
    
    

    