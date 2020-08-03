from igraph import *
import pyUtils

# # Remote debugger setup (KEEP IT COMMENTED UNLESS DEBUGGING!)
# import sys
# sys.path.append("C:\\Users\\ityuryukanov\\.p2\\pool\\plugins\\org.python.pydev_4.3.0.201508182223\\pysrc")
# import pydevd  
# pydevd.settrace()                                                                               

def mat2igraph(edges, vs_lbl=None, es_grp=None, vs_grp=None):
    """ Create an UNDIRECTED graph according to the input specification, making use of the
        rich capabilities of the python-igraph library.
        
        edges is a list of (i, j, w) tuples: i, j are the vertex NAMES (g.vs['name']) of
              an edge (i, j are from 0 to num_nodes - 1) and w is the edge weight
        vs_lbl is a list of node labels: g.vs[i]['label'] == vs_lbl(g.vs[i]['name']), and
              g.vs[i]['name'] comes from the list edges which specifies the graph
        es_grp is a list of edges which should be colored in a distinctive color (red): 
              the two-tuples in the list represent g.vs['name'] values of respective nodes
        vs_grp is a list of tuples of vertex NAMES (g.vs['name']) which belong to the same
              group, or the list of colors of each vertex (color#1 for vertex '1' and so on)              
        
        Requires: python igraph v.0.7.0 and higher
        
        Author: Ilya Tyuryukanov
        Date: 30 September 2015
        Last revision: 2 October 2015 """   
    # Create weighted graph from 'edges' (i, j, w) - list 
    g = Graph.TupleList(edges, directed=False, weights=True)
    num_nodes = len(g.vs)
    num_edges = len(g.es)
    
    # Assign node labels (e.g. bus numbers from the initital model) for plotting
    if vs_lbl is not None:
        g.vs['label'] = [str(vs_lbl[g.vs[i]['name']]) for i in range(0,len(vs_lbl))]
    
    # Assign edge width (for drawing) proportional to edge weight
    wght = [w for _, _, w in edges]
    g.es['label'] = ['%.2f'%w for w in wght]
    g.es['width'] = [1 for _, _, w in edges]
    g.es['weight'] = [1 for _, _, w in edges]
    # wght_min = min(wght) 
    # wght_range = max(wght) - wght_min
    # ES_MIN = 1  # min line width
    # ES_MAX = 1.1  # max line width
    # es_range = ES_MAX - ES_MIN
    # es_wdth = [(w - wght_min) / wght_range * es_range + ES_MIN for w in wght]  # stackoverflow.com  
    # g.es['width'] = es_wdth
    
    # Change color for edges in 'es_grp' list to red. Here order of vertices in an edge 
    # may be important: ROBUST check of BOTH orders due to UNDIRECTED graph assumption
    EDGE_COL = '#0A0A0A'  # default edge color
    RED = '#FF0000' 
    es_col = list()
    if es_grp is not None:
        for e in g.es:
            curr_edge1 = (g.vs[e.source]['name'], g.vs[e.target]['name'])
            curr_edge2 = (g.vs[e.target]['name'], g.vs[e.source]['name'])
            if (curr_edge1 in es_grp) or (curr_edge2 in es_grp):
                es_col.append(RED)
            else:
                es_col.append(EDGE_COL)
    else:
        es_col = [EDGE_COL]*num_edges
    g.es['color'] = es_col            
    
    # Set node colors
    NODE_COL = '#AAAAAA'  # default node color rgb(170,170,170)
    vs_col = [NODE_COL]*num_nodes   
    if vs_grp is not None:    
        tpl_len = [len(tpl) for tpl in vs_grp]   
        log_len = [True if num==3 else False for num in tpl_len]
        tpl_val = [True if item<1.00001 and item>-0.00001 else False for item in list(sum(vs_grp, ()))] 
        if all(log_len) and all(tpl_val):  # If input is rgb-color tuples for each node
            vs_grp = rgb_to_hex(vs_grp)  # convert rgb color attributes to hex
            for i in range(0, num_nodes):
                vs_col[i] = vs_grp[g.vs[i]['name']]
        else:  # If input is a list of tuples of node groups (without colors)
            n_grp = len(vs_grp)  # number of node groups to be colored
            grp_col = get_n_col(n_grp)
            for i in range(0, num_nodes):
                for j in range(0, n_grp):
                    if g.vs[i]['name'] in vs_grp[j]:
                        vs_col[i] = grp_col[j]                    
    g.vs['color'] = vs_col
    g.vs['fillcolor'] = vs_col
    return g
    
def plot_igraph(g, file, bbox, layout, margin, v_sz = 10, v_lbl_sz = 9, v_shape = 'circle', 
                v_lbl = True, e_lbl = True):
    """ Saves an igraph object with predefined attributes using a chosen layout and 
        margin, with or without node/edge labels to an svg file """  
    # Enforce choice between default labels or no labels
    if v_lbl:
        v_lbl = g.vs['label']
    else:
        v_lbl = None
        
    if e_lbl:
        e_lbl = g.es['label']
    else:
        e_lbl = None    
        
    plot(g, target = file+'.svg', bbox = bbox, layout = layout, margin = margin, 
         vertex_label = v_lbl, edge_label = e_lbl, vertex_size = v_sz, 
         vertex_label_size = v_lbl_sz, vertex_shape = v_shape)
         
         
def delete_edges(g, e_list):
    """ Delete edges with g.vs['name'] equal to two-tuples from e_list """
    if not isinstance(e_list, list) or len(e_list) == 0 or len(e_list[0]) != 2  \
      or not isinstance(e_list[0], tuple) or not isinstance(e_list[0][0], int): 
        return g # later throw an error here!!!
        
    if not isinstance(g, Graph):
        return g # later throw an error here!!!
    
    e_idx = list()
    for e in g.es:
        curr_edge1 = (g.vs[e.source]['name'], g.vs[e.target]['name'])
        curr_edge2 = (g.vs[e.target]['name'], g.vs[e.source]['name'])
        if (curr_edge1 in e_list) or (curr_edge2 in e_list):
            e_idx.append(e.index)

    g.delete_edges(e_idx)      
    return g    
        

def igraph2graphviz(g, file):
    """" Save an igraph object with predefined attributes to the GraphViz format. 
         A recommended way to plot large graphs is by:
         sfdp -Gsize=67! -Goverlap=prism -Tpng root.gv > root.png """
    
    g.write_dot(file+'.dot')
    f = open(file+'.dot', "r")
    contents = f.readlines()
    f.close()    
    contents = pyUtils.strsci2strfix(contents)
    contents.insert(2, 'node [ style = filled ];\n')
    f = open(file+'.dot', "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()


def get_n_col(num_colors=5):
    """ Gen num_colors distinct colors on the HSV scale with s = 0.4, v = 0.99 (bright)
        Range mapping is only valid for colors problem (i.e. map from [0,1] to [0,255])
        because 255 is a suffiently large number to hide some intrinsic inaccuracy
        From http://stackoverflow.com/ """
    import colorsys
    
    hsv_tuples = [(x*1.0/num_colors, 0.4, 0.99) for x in range(num_colors)]
    rgb_tuples = list(map(lambda x: colorsys.hsv_to_rgb(*x), hsv_tuples))
    hex_colors = rgb_to_hex(rgb_tuples)
    return hex_colors 
    
    
def rgb_to_hex(rgb_tuples):
    hex_colors = ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) for (r,g,b) in rgb_tuples]          
    return hex_colors     
    