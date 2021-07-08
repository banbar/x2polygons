"""
This module contains the functions to calculate the distance between two matching polygons (e.g. building footprints coming from two different datasets). 

Use cases:
    The functionality presented in this module could be used to compare reference building datasets with 
    a VGI dataset such as OSM.
    
Distance Functions:
    - ``Chamfer Distance`` 
    - ``Hausdorff Distance``
    - ``PoLis Distance`` `[1] <https://ieeexplore.ieee.org/document/6849454>`_.
    - ``Turn Function Distance`` `[2] <https://ieeexplore.ieee.org/document/75509>`_.


"""

import math
import copy
from shapely.geometry import Polygon, Point
import geopandas as gp

# When packaging & developing:
#from . import plot as plt

# When creating the documentation
import plot as plt


def polygon_vertices(polygon):
    '''
    Converts an input polygon into its nodes.
    
    Args:
        - **polygon** (*polygon*): A polygon object

    Returns:
        - **point_list** (*point []*):  The ordered set of nodes of the input polygon
    '''
    
    pointsx = list(list(polygon.exterior.coords.xy)[0]) # X-axis coordinates
    pointsy = list(list(polygon.exterior.coords.xy)[1]) # Y-axis coordinates
    point_list = []
    for q in range(len(pointsx)):
        point_list.append([pointsx[q],pointsy[q]])
    return point_list

class point:
    """
    A class to represent a point.
    
    Args:
        - **x** (*float*): The x coordinate of a point
        - **y** (*float*): The y coordinate of a point
    
    Methods:
        - **distance_to_point** (*px*): Returns the distance to the point *px*. 
    
    """
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def distance_to_point(self, px):
        return ( (self.x-px.x)**2 + (self.y-px.y)**2 )
        
    
class line_vector:
    """
    A class to represent a line vector composed of two points. 
    
    Args:
        - **p1** (*point*): Start point of the vector
        - **p2** (*point*): End point of the vector
    
    Methods:
        - **point_on_where** (*px*): Identifies the orientation of the point px with respect to the vector. Returns 'LEFT', 'RIGHT' or 'COLINEAR'.
        - ** length**: Returns the length of the vector
        - **angle_to_vector** (*vx*): Returns the degree between the self vector and the query vector vx.
    
   
    """
    def __init__(self, p1, p2):
        # the vector is from point 1 (p1) to p2
        self.p1 = p1
        self.p2 = p2
    
    def point_on_where(self, px): 
       
        # First vector is from p1 -> p2
        v1 = point(0, 0)
        v2 = point(0, 0)
        
        v1.x = self.p2.x - self.p1.x 
        v1.y = self.p2.y - self.p1.y
        
        v2.x = px.x - self.p2.x
        v2.y = px.y - self.p2.y
        
        # Find the determinant - Right hand rule
        result = v1.x*v2.y - v1.y*v2.x
        
        if(result > 0):
            return 'LEFT' # yes, px is on the LEFT hand-side of the line
        elif (result < 0):
            return 'RIGHT'
        else: # result = 0 - point is colinear with the line
            return 'COLINEAR'
        
        
    def length(self):
        return ( math.sqrt( (self.p2.x-self.p1.x)**2 + (self.p2.y-self.p1.y)**2 ) )
    
    def angle_to_vector(self, vx):

        # https://onlinemschool.com/math/library/vector/angl/#:~:text=Definition.,Basic%20relation.
        cos_alpha = ( (self.p2.x-self.p1.x)*(vx.p2.x - vx.p1.x) + (self.p2.y-self.p1.y)*(vx.p2.y - vx.p1.y)) / (self.length()*vx.length()) 
        
        radian = math.acos(cos_alpha)
        degree = math.degrees(radian) # convert radian to degrees
    
        return degree

def chamfer_distance(polygon_a, polygon_b, **kwargs):
    '''
    Identifies the Chamfer distance between two input polygons. The distance is calculated from *polygon_a* to *polygon_b* (a->b).
    
    Args:
        - **polygon_a** (*polygon*): First polygon
        - **polygon_b** (*polygon*): Second polygon
        - **kwargs**:
            - symmetrise: How to symmetrise the distance measure as there would be two distances (i.e. a->b, b->a). Options are: *'min'*, *'max'*, *'average'*. The *average* option is weighted by the number of nodes of each polygon as described `here <https://ieeexplore.ieee.org/document/6849454>`_.
            
    Returns:
        - **distance** (*float*): Chamfer distance between the polygons
    '''
    
    c_a_b = 0 # init the directed Chamfer Distance between polygon A and B
    c_b_a = 0 # init the directed Chamfer Distance between polygon B and A
    
    vertices_a = polygon_vertices(polygon_a)
    vertices_b = polygon_vertices(polygon_b)
    
    for i in range(len(vertices_a)-1): # from each corner of the polygon 1
        minimum_distance = 1000.0 # Minimum distance set as initial.
        for j in range(len(vertices_b)): # to each corner of the polygon 2
            distance = ((vertices_a[i][0] - vertices_b[j][0])**2+(vertices_a[i][1] - vertices_b[j][1])**2)**0.5 # The distance between corners is calculated
            if minimum_distance > distance: # If the calculated distance is greater than the minimum distance
                minimum_distance = distance # minimum distance is calculated distance
        c_a_b += minimum_distance # Add minimum distance to total distance
     
    
    for k in range(len(vertices_b)-1): # from each corner of the polygon 2
        minimum_distance = 1000.0 # Minimum distance set as initial.
        for l in range(len(vertices_a)): # to each corner of the polygon 2
            distance = ((vertices_a[l][0] - vertices_b[k][0])**2+(vertices_a[l][1] - vertices_b[k][1])**2)**0.5 # The distance between corners is calculated
            if minimum_distance > distance: # If the calculated distance is greater than the minimum distance
                minimum_distance = distance # minimum distance is calculated distance
        c_b_a += minimum_distance # Add minimum distance to total distance
    
    # Default: c_a_b
    if('symmetrise' not in kwargs):
        return c_a_b
    elif(kwargs['symmetrise'] == 'average'):
        return ( (c_a_b / (2* (len(vertices_a)-1)) ) + (c_b_a / (2* (len(vertices_b)-1)) ) )
    elif(kwargs['symmetrise'] == 'min'):
        return min(c_a_b, c_b_a)
    elif(kwargs['symmetrise'] == 'max'):
        return max(c_a_b, c_b_a)
    


    
def hausdorff_distance(polygon_a, polygon_b, **kwargs):
    '''
    Identifies the Hausdorff distance between two input polygons. The distance is calculated from *polygon_a* to *polygon_b* (a->b).
    
    Args:
        - **polygon_a** (*polygon*): First polygon
        - **polygon_b** (*polygon*): Second polygon
        - **kwargs**:
            - symmetrise: How to symmetrise the distance measure as there would be two distances (i.e. a->b, b->a). Options are: *'min'*, *'max'*, *'average'*. 
            
    Returns:
        - **distance** (*float*): Hausdorf distance between the polygons
    '''
    
    distance_between_vertices = []
    
    vertices_a = polygon_vertices(polygon_a)
    vertices_b = polygon_vertices(polygon_b)
    
    
    for i in range(len(vertices_a)): # from each corner of the polygon 1
        minimum_distance = 1000.0 # Minimum distance set as initial.
        for j in range(len(vertices_b)): # to each corner of the polygon 2
            distance = ((vertices_a[i][0] - vertices_b[j][0])**2+(vertices_a[i][1] - vertices_b[j][1])**2)**0.5 # The distance between corners is calculated
            if minimum_distance > distance: # If the calculated distance is greater than the minimum distance
                minimum_distance = distance # minimum distance is calculated distance
        distance_between_vertices.append(minimum_distance) # the minimum distance is added to the list
    h_a_b = max(distance_between_vertices) # The greatest value between the smallest distances becomes the Hausdorff distance
    
    distance_between_vertices = []
    
    for k in range(len(vertices_b)): # from each corner of the polygon a
        minimum_distance = 1000.0 # Minimum distance set as initial.
        for l in range(len(vertices_a)): # to each corner of the polygon a
            distance = ((vertices_a[l][0] - vertices_b[k][0])**2+(vertices_a[l][1] - vertices_b[k][1])**2)**0.5 # The distance between corners is calculated
            if minimum_distance > distance: #  If the calculated distance is greater than the minimum distance
                minimum_distance = distance # minimum distance is calculated distance
        distance_between_vertices.append(minimum_distance) # the minimum distance is added to the list
    h_b_a = max(distance_between_vertices) # The greatest value between the smallest distances becomes the Hausdorff distance
    
    # default options:
        # directed = False
        # symmetrize = max
    if('symmetrise' not in kwargs):
        return h_a_b
    elif(kwargs['symmetrise'] == 'min'):
        return min(h_a_b, h_b_a) # Number of nodes may DIFFER -  OSM may have many nodes on the line - shapes are very similar but Hausdorrf distance is large - if we were optimist, then the H distance could have been zero
    elif(kwargs['symmetrise'] == 'max'):
        return max(h_a_b, h_b_a)
    elif(kwargs['symmetrise'] == 'average'):
        return (h_a_b + h_b_a)/2

def polis_distance(polygon_a, polygon_b, **kwargs):
    '''
    Identifies the PoLis distance between two input polygons. The distance is calculated from *polygon_a* to *polygon_b* (a->b).
    
    Args:
        - **polygon_a** (*polygon*): First polygon
        - **polygon_b** (*polygon*): Second polygon
        - **kwargs**:
            - symmetrise: How to symmetrise the distance measure as there would be two distances (i.e. a->b, b->a). Options are: *'min'*, *'max'*, *'average'*. 
            
    Returns:
        - **distance** (*float*): PoLis distance between the polygons
    '''
    geoSeriesA = gp.GeoSeries(polygon_a)
    geoSeriesB = gp.GeoSeries(polygon_b)
    # We can hold a VISITED polygon list - we can skip those to improve the run-time

    # Convert the polygon into a geoseries onject
    #print(geoseries_B)

    # vertices of A  -> to -> polygon B
    # For all the vertices of A:
    num_vertices_a = len(polygon_a.exterior.coords) - 1 # obtain the number of vertices - first & last vertices coincide
    distance_to_edges_a_b = 0

    for i in range(num_vertices_a):
        # Convert each vertice to a GeoSeries point object:
        tmp_vertex = gp.GeoSeries(Point(polygon_a.exterior.coords[i]))

        # Calculate the distance between the point and the polygon B
        dist = geoSeriesB.boundary.distance(tmp_vertex) # distance function does not work properly - if a vertex is inside the other polygon, it will return zero regardless its position
        distance_to_edges_a_b += dist[0]
    
    polis_a_b = distance_to_edges_a_b / num_vertices_a

    # vertices of B  -> to -> polygon A
    # For all the vertices of A:
    num_vertices_b = len(polygon_b.exterior.coords) - 1  # obtain the number of vertices
    distance_to_edges_b_a = 0

    for j in range(num_vertices_b):
        # Convert each vertice to a GeoSeries point object:
        tmp_vertex = gp.GeoSeries(Point(polygon_b.exterior.coords[j]))

        # Calculate the distance between the point and the polygon A
        dist = geoSeriesA.boundary.distance(tmp_vertex) # distance function does not work properly - if a vertex is inside the other polygon, it will return zero regardless its position
        distance_to_edges_b_a += dist[0]
    
    polis_b_a = distance_to_edges_b_a / num_vertices_b
    
    vertices_a = polygon_vertices(polygon_a)
    vertices_b = polygon_vertices(polygon_b)

    # Calculate PoLiS
    # Default: polis_a_b (directed)
    if('symmetrise' not in kwargs):
        return polis_a_b
    elif(kwargs['symmetrise'] == 'average'):
        return ( (polis_a_b / 2) + (polis_b_a / 2) )
    elif(kwargs['symmetrise'] == 'max'):
        return max(polis_a_b, polis_b_a)
    elif(kwargs['symmetrise'] == 'min'):
        return min(polis_a_b, polis_b_a)

    



def turn_function(polygon, **kwargs):
    '''
    Identifies the turn function of an input polygon. 
    
    Args:
        - **polygon** (*polygon*): The input polygon
    
    Returns:
        - **dictionary** with the following attributes
            - **angles** (*float []*): turn angles 
            - **lengths** (*float []*): normalised lengths
            - **direction** (*char []*): each turn direction - values of the list are Left (*l*), Right (*R*) or Colinear (*-*)
            - **digitisation_direction** (*str*): Counter ClockWise (*CCW*) or ClockWise (*CW*)
        - **kwargs**:
            - **ccw**: Enforce a ccw turn (*True* or *False* (default))
            - **plot**: Plot the turn function of the polygon (*True* or *False* (default))
    '''   
   
    tmp_points = polygon_vertices(polygon)
    points = []
    # save it as a point
    for i in range(len(tmp_points)):
        points.append(point(tmp_points[i][0], tmp_points[i][1] ))
    
    v_init = line_vector(points[0], points[1])
    
    # For all the remaining vectors 
    # 1. Identify the change in direction - Counter-Clock Wise (CCW, Left) or Clock Wise (CW, Right)
    # 2. Find the angle between the vectors
    # 3. Total is + by CCW changes, - by CW changes
    
    total_length = 0
    turn = {}
    turn['angles'] = []
    turn['lengths'] = []
    turn['direction'] = []
    turn['digitisation_direction'] = 'CCW' # Assume CCW: +360 
        
    for i in range(len(points)-2):
        vx = line_vector(points[i], points[i+1])
        vy = line_vector(points[i+1], points[i+2])
        
        # Start adding the length of segments from the 2. segment
        # Then we will add the first segment's length
        if(vx.point_on_where(points[i+2]) == 'LEFT'):
            turn["angles"].append(vx.angle_to_vector(vy))
            turn["direction"].append('L') 
        elif(vx.point_on_where(points[i+2]) == 'RIGHT'):
            turn["angles"].append(-vx.angle_to_vector(vy))
            turn["direction"].append('R') 
        else: # COLINEARITY - do not make a turn !!!
            turn["angles"].append(0)            
            turn["direction"].append('-') 
            
        turn["lengths"].append(vy.length())
        
    # Add the angle between the last and first segments
    if(vy.point_on_where(v_init.p2) == 'LEFT'):
        turn["angles"].append(vy.angle_to_vector(v_init))
        turn["direction"].append('L') 
    elif (vy.point_on_where(v_init.p2) == 'RIGHT'):
        turn["angles"].append(-vy.angle_to_vector(v_init))
        turn["direction"].append('R') 
    else: # COLINEAR
        turn["angles"].append(0)            
        turn["direction"].append('-') 
    
    turn["lengths"].append(v_init.length())
        
    # Normalise the lengths
    total_length = sum(turn["lengths"])
    for i in range(len(turn["lengths"])):
        turn["lengths"][i] = turn["lengths"][i] / total_length
    
    # Post-Process the COLINEARITY
    indices_to_delete = []
    for i in range(len(turn["lengths"])-1, 0, -1):
        if(turn["direction"][i] == '-'):
            turn["lengths"][i-1] += turn["lengths"][i]
            indices_to_delete.append(i)
    
    # Handle the first change
    found = 0
    if(turn["direction"][0] == '-'):
        indices_to_delete.append(0)
        while(not found):
            for i in range(len(turn["lengths"])-1, 0, -1):
                if(turn["direction"][i] != '-'):
                    found = 1
                    break
    if(found):
        turn["lengths"][i] += turn["lengths"][0]
                    
    # Reorganise the turn function 
    sorted_indecies_to_delete = sorted(indices_to_delete, reverse=True)
    
    if(len(indices_to_delete) > 0):
        for i in (indices_to_delete):
            del turn["direction"][i]
            del turn["lengths"][i]
            del turn["angles"][i]
   
    
    # Output a turn function by moving towards the reverse order
    #!!!!!!!!!!!!!!!    CCW = TRUE !!!!!!!!!!!!!!!
    
    if( round(sum(turn["angles"])) == -360):
        turn['digitisation_direction'] = 'CW' 
        if ('ccw' in kwargs):
            v_init = line_vector(points[-1], points[-2])
            
            turn['angles'] = []
            turn['lengths'] = []
            turn['direction'] = []
    
                    
            total_length = 0
               
            for i in range(len(points)-2):
                vx = line_vector(points[-i-1], points[-i-2])
                vy = line_vector(points[-i-2], points[-i-3])
                
                # Determine the turn
                if(vx.point_on_where(points[-i-3]) == 'LEFT'):
                    turn["angles"].append(vx.angle_to_vector(vy))
                    turn["direction"].append('L') 
                elif(vx.point_on_where(points[-i-3]) == 'RIGHT'):
                    turn["angles"].append(-vx.angle_to_vector(vy))
                    turn["direction"].append('R') 
                else: # COLINEARITY - do not make a turn !!!
                    turn["angles"].append(0)     
                    turn["direction"].append('-') 
                
                turn["lengths"].append(vy.length())
            
            # Add the angle between the last and first segments
            if(vy.point_on_where(v_init.p2) == 'LEFT'):
                turn["angles"].append(vy.angle_to_vector(v_init))
                turn["direction"].append('L') 
            elif (vy.point_on_where(v_init.p2) == 'RIGHT'):
                turn["angles"].append(-vy.angle_to_vector(v_init))
                turn["direction"].append('R') 
            else: # COLINEAR
                turn["angles"].append(0)            
                turn["direction"].append('-') 
                
            turn["lengths"].append(v_init.length())
                    
            # Normalise the lengths
            total_length = sum(turn["lengths"])
            for i in range(len(turn["lengths"])):
                turn["lengths"][i] = turn["lengths"][i] / total_length


            # Post-Process the COLINEARITY
            indices_to_delete = []
            for i in range(len(turn["lengths"])-1, 0, -1):
                if(turn["direction"][i] == '-'):
                    turn["lengths"][i-1] += turn["lengths"][i]
                    indices_to_delete.append(i)
            
            # Handle the first change
            found = 0
            if(turn["direction"][0] == '-'):
                indices_to_delete.append(0)
                while(not found):
                    for i in range(len(turn["lengths"])-1, 0, -1):
                        if(turn["direction"][i] != '-'):
                            found = 1
                            break
            if(found):
                turn["lengths"][i] += turn["lengths"][0]
                    

            sorted_indecies_to_delete = sorted(indices_to_delete, reverse=True)
            
            if(len(indices_to_delete) > 0):
                for i in (indices_to_delete):
                    del turn["direction"][i]
                    del turn["lengths"][i]
                    del turn["angles"][i]
    
    if ('plot' in kwargs):
        plt.plot_turn_function(turn)
        
    return turn



def distance_between_turn_functions(a_turn, b_turn):
    '''
    Calculates the distance between two turn functions. 
    
    Args:
        - **a_turn** (*dict*): First polygon's turn function as provided by the *turn_function(polygon_a)*.
        - **b_turn** (*dict*): Second polygon's turn function as provided by the *turn_function(polygon_b)*.
    
    Returns:
        - **dict**: A dictionary containing the following attributes:
            - **a** (*dict*): aligned turn function of polygon a
            - **b** (*dict*): aligned turn function of polygon b
            - **distance** (*float*): turn function distance between the polygons
            - **piece_wise_a** (*float []*): piece wise lengths of polygon a
            - **piece_wise_b** (*float []*): piece wise lengths of polygon b
            - **combined_piece_wise** (*float []*): combined piece wise lengths of two polygons
    Notes:
        - Alignment step is carried out to make sure that the start point of the polygons do not make a difference.
        
        

    '''   
    # Notes:
        # Number of edges may differ between two polygons
    #1. Align the turn angles so that the distance is minimum
    # - We may start from different nodes, the node numbers may be different etc.
    # pop & append
    distances = []
    min_distance = 10**10
    
    # Normalise the angles
    # Round the float - we may see error in floating point arithmetic
    total_angle = sum(a_turn["angles"]) # must be 360
    for i in range(len(a_turn["angles"])):
        a_turn["angles"][i] = round((a_turn["angles"][i] / total_angle), 3)
    
    for i in range(len(b_turn["angles"])):
        b_turn["angles"][i] = round((b_turn["angles"][i] / total_angle), 3)
    
    # Round the lengths
    for i in range(len(a_turn["lengths"])):
        a_turn["lengths"][i] = round(a_turn["lengths"][i], 3)
    
    for i in range(len(b_turn["lengths"])):
        b_turn["lengths"][i] = round(b_turn["lengths"][i], 3)
        
    
    
    
    
    piece_wise_lengths = {} 
    piece_wise_lengths["b"] = [0]*(len(b_turn["lengths"])+1)
    for i in range(1, len(b_turn["lengths"])+1):
        piece_wise_lengths["b"][i] = piece_wise_lengths["b"][i-1] + b_turn["lengths"][i-1]
        
    # SHIFT Polygon A: -----------------------------------------    
    for shift in range(len(a_turn["lengths"])):     
        piece_wise_lengths["a"] = [0]*(len(a_turn["lengths"])+1)
        for i in range(1, len(a_turn["lengths"])+1):
            piece_wise_lengths["a"][i] = piece_wise_lengths["a"][i-1] + a_turn["lengths"][i-1]
        
    
        combined_piece_wise_lengths = [] + piece_wise_lengths["a"] + piece_wise_lengths["b"]
        combined_piece_wise_lengths.sort()
        # remove the first and last elements - i.e. 0 & 1
        combined_piece_wise_lengths.pop(0)
        combined_piece_wise_lengths.pop(len(combined_piece_wise_lengths)-1)
        
        # Round the combined_piece_wise_lengths - strange why we see still not rounded values
        for i in range(len(combined_piece_wise_lengths)):
            combined_piece_wise_lengths[i] = round(combined_piece_wise_lengths[i], 3)
        
        distance = 0
        index_a = 1
        index_b = 1
        
        for i in range(1, len(combined_piece_wise_lengths)):
            for j in range(1, len(piece_wise_lengths["a"])):
                if(piece_wise_lengths["a"][j-1] < combined_piece_wise_lengths[i] <= piece_wise_lengths["a"][j] ):
                    # use j-1
                    index_a = j-1
                    break
            
            for k in range(1, len(piece_wise_lengths["b"])):
                if(piece_wise_lengths["b"][k-1] < combined_piece_wise_lengths[i] <= piece_wise_lengths["b"][k] ):
                    # use j-1
                    index_b = k-1
                    break
            
            distance += abs(a_turn["angles"][index_a] - b_turn["angles"][index_b])
            
        if(distance < min_distance): # save a snapshot
            min_distance_snapshot = {}
            min_distance = distance
            min_distance_snapshot["a"] = copy.deepcopy(a_turn)
            min_distance_snapshot["b"] = copy.deepcopy(b_turn)
            min_distance_snapshot["distance"] = min_distance
            min_distance_snapshot["piece_wise_a"] = copy.deepcopy(piece_wise_lengths["a"])
            min_distance_snapshot["piece_wise_b"] = copy.deepcopy(piece_wise_lengths["b"])
            min_distance_snapshot["combined_piece_wise"] = []+combined_piece_wise_lengths
        
        distances.append(distance)
        
        # SHIFT Polygon A: 
        angle = a_turn["angles"].pop(0)
        length = a_turn["lengths"].pop(0)
        a_turn["angles"].append(angle)
        a_turn["lengths"].append(length)
    
    # SHIFT Polygon B: -----------------------------------------  
    
    piece_wise_lengths = {} 
    piece_wise_lengths["a"] = [0]*(len(a_turn["lengths"])+1)
    
    for i in range(1, len(a_turn["lengths"])+1):
        piece_wise_lengths["a"][i] = piece_wise_lengths["a"][i-1] + a_turn["lengths"][i-1]
        
    for shift in range(len(b_turn["lengths"])):     
        piece_wise_lengths["b"] = [0]*(len(b_turn["lengths"])+1)
        for i in range(1, len(b_turn["lengths"])+1):
            piece_wise_lengths["b"][i] = piece_wise_lengths["b"][i-1] + b_turn["lengths"][i-1]
        
    
        combined_piece_wise_lengths = [] + piece_wise_lengths["a"] + piece_wise_lengths["b"]
        combined_piece_wise_lengths.sort()
        # remove the first and last elements - i.e. 0 & 1
        combined_piece_wise_lengths.pop(0)
        combined_piece_wise_lengths.pop(len(combined_piece_wise_lengths)-1)
        
        distance = 0
        index_a = 1
        index_b = 1
        
        for i in range(1, len(combined_piece_wise_lengths)):
            for j in range(1, len(piece_wise_lengths["a"])):
                if(piece_wise_lengths["a"][j-1] < combined_piece_wise_lengths[i] <= piece_wise_lengths["a"][j] ):
                    # use j-1
                    index_a = j-1
                    break
            
            for k in range(1, len(piece_wise_lengths["b"])):
                if(piece_wise_lengths["b"][k-1] < combined_piece_wise_lengths[i] <= piece_wise_lengths["b"][k] ):
                    # use j-1
                    index_b = k-1
                    break
            
            distance += abs(a_turn["angles"][index_a] - b_turn["angles"][index_b])
            
        if(distance < min_distance): # save a snapshot
            min_distance_snapshot = {}
            min_distance = distance
            min_distance_snapshot["a"] = copy.deepcopy(a_turn)
            min_distance_snapshot["b"] = copy.deepcopy(b_turn)
            min_distance_snapshot["distance"] = min_distance
            min_distance_snapshot["piece_wise_a"] = copy.deepcopy(piece_wise_lengths["a"])
            min_distance_snapshot["piece_wise_b"] = copy.deepcopy(piece_wise_lengths["b"])
            min_distance_snapshot["combined_piece_wise"] = []+combined_piece_wise_lengths

            
        distances.append(distance)
        
        # SHIFT Polygon B: 
        angle = b_turn["angles"].pop(0)
        length = b_turn["lengths"].pop(0)
        b_turn["angles"].append(angle)
        b_turn["lengths"].append(length)
    
    # Calculate the total distance of the min_distance_snapshot
    index_a = 1
    index_b = 1
    total_distance = 0
    # remove the possibly repeating change points
    min_distance_snapshot["combined_piece_wise"] = sorted(set(min_distance_snapshot["combined_piece_wise"]))
    for i in range(1, (len(min_distance_snapshot["combined_piece_wise"]))):
        for j in range(index_a, len(min_distance_snapshot["piece_wise_a"])):
            if(min_distance_snapshot["piece_wise_a"][j-1] < min_distance_snapshot["combined_piece_wise"][i] <= min_distance_snapshot["piece_wise_a"][j] ):
                # use j-1
                index_a = j-1
                break
            
        for k in range(index_b, len(min_distance_snapshot["piece_wise_b"])):
            if(min_distance_snapshot["piece_wise_b"][k-1] < min_distance_snapshot["combined_piece_wise"][i] <= min_distance_snapshot["piece_wise_b"][k] ):
                # use j-1
                index_b = k-1
                break
        
        total_distance += abs(min_distance_snapshot["a"]["angles"][index_a] - min_distance_snapshot["b"]["angles"][index_b]) + abs(min_distance_snapshot["a"]["lengths"][index_a] - min_distance_snapshot["b"]["lengths"][index_b])
 
    
    min_distance_snapshot["total_distance"] = total_distance
        
        
    
    return min_distance_snapshot

def turn_function_distance(p1, p2, **kwargs):
    # Obtain the turn functions
    if('ccw' in kwargs):
        p1_turn = turn_function(p1, ccw = kwargs["ccw"])
        p2_turn = turn_function(p2, ccw = kwargs["ccw"])
    else:
        p1_turn = turn_function(p1)
        p2_turn = turn_function(p2)
    
    return distance_between_turn_functions(p1_turn, p2_turn)["total_distance"]

    
