"""Plot functions. At the moment only Turn Function"""

import matplotlib.pyplot as plt

def plot_turn_function(turn):  
    '''
    Plots the cumulative length turn function. 
    
    Args:
        turn (dict): The turn dictionary of a polygon.
    
    Returns:
        plot
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