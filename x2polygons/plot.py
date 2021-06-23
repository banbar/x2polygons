""" This module contains all plot related functionality.
"""

import matplotlib.pyplot as plt

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
                 [poly_a.exterior.coords[i][1], poly_a.exterior.coords[i+1][1]], 'b')
        
    # Plot the nodes
    for i in range(len(poly_a.exterior.coords)-1):
        ax.plot(poly_a.exterior.coords[i][0], 
                 poly_a.exterior.coords[i][1], 'bo',
                 markersize = 12,
                 fillstyle = 'none')
    
    # Second polygon
    # Plot the edges 
    for i in range(len(poly_b.exterior.coords)-1):
        ax.plot([poly_b.exterior.coords[i][0], poly_b.exterior.coords[i+1][0]], 
                 [poly_b.exterior.coords[i][1], poly_b.exterior.coords[i+1][1]], 'r',
                 linestyle='dashed')
        
    # Plot the nodes
    for i in range(len(poly_b.exterior.coords)-1):
        ax.plot(poly_b.exterior.coords[i][0], 
                 poly_b.exterior.coords[i][1], 'rs',
                 markersize = 12,
                 fillstyle = 'none')
    
    # Remove the axes
    fig.patch.set_visible(False)
    ax.axis('off')
    
    if('save_plot_as' in kwargs):
        fig.savefig(kwargs["save_plot_as"], format = 'svg', dpi=300)
    

    
    
    