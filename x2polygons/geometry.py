# -*- coding: utf-8 -*-
"""
This module has the geometry related functions. The outcomes of the functions used in this module can be used to calculate accuracy, F-1 score etc. Therefore, it only aoperates on 1-1 matching building footprints. 

"""

import math
import copy
from shapely.geometry import Polygon, Point
import geopandas as gp

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
        - **point_on_where** (self, *px*): Identifies the orientation of the point *px* with respect to the vector. Returns 'LEFT', 'RIGHT' or 'COLINEAR'.
        - ** length**(self): Returns the length of the vector
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
    

def x2_areas(polygon_test, polygon_ref):
    '''
    Identifies the intersecting area of two input (geopandas) polygons (True Positive, TP), False-Positive (FP, area belonging to test (e.g. OSM) but not reference) and False-Negative (FN, area belonging to reference polygon but not test). The function operates only on 1-1 matching polygons.
    
    Args:
        - **polygon_test** (*polygon*): Test (e.g. OSM) polygon
        - **polygon_ref** (*polygon*): Reference polygon 
            
    Returns:
        - **result** (*dict*): Dictionary reporting the three metrics: TP, FP, FN
    
    .. image:: https://raw.githubusercontent.com/banbar/x2polygons/ae626c5ee9042b710b913099b04693b751274b28/img/area_overlap.png
       :width: 300px
       :height: 300px
    
    '''
    
    result = {}
    
    result["TP"] = polygon_ref.intersection(polygon_test).area
    result["FP"] = polygon_test.area - result["TP"]
    result["FN"] = polygon_ref.area - result["TP"]
    
    
    return result

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

def polygon_perimeter(polygon):
    '''
    Returns the perimeter of a polygon. 
    
    Args:
        - **polygon** (*polygon*): A polygon object

    Returns:
        - **perimeter** (*float*):  Perimeter of the input polygon. 
    '''
    
    pointsx = list(list(polygon.exterior.coords.xy)[0]) # X-axis coordinates
    pointsy = list(list(polygon.exterior.coords.xy)[1]) # Y-axis coordinates
    
    perimeter = 0
    for i in range(len(pointsx)-1):
        p1 = point(pointsx[i], pointsy[i])
        p2 = point(pointsx[i+1], pointsy[i+1])
        
        # Convert to an edge 
        edge = line_vector(p1, p2)
        
        perimeter += edge.length()
        
    return perimeter

