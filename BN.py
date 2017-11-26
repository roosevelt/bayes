# Bayesian network class 
from Node import Node
from Edge import Edge
import pydot

class BayesianNetwork:
    
    def __init__(self, nodes):
        self.nodes =  nodes
        self.edges = self.__buildEdges(self.nodes)
        
    def __buildEdges(self, nodes):
        edge_list = []
        for node in nodes:
            parents = node.CPD.condition_variables
            for parent in parents:
                new_edge = Edge(parent, node)
                edge_list.append(new_edge)
        return edge_list
    
    def __str__(self):
        n_repetitions = 20
        buffer=""
        buffer+="*"*n_repetitions + " Bayesian Network Description " + "*"*n_repetitions + "\n"
        buffer+="Nodes: " + ", ".join([node.name for node in self.nodes]) + "\n"
        buffer+="Number of Edges: " + str(len(self.edges)) + "\n\n"
        buffer+="Conditional Probability Tables (all nodes): \n\n"
        for node in self.nodes:
            buffer+=node.__str__()+ "\n"
        return buffer
            
    def exportGraph(self, path="bayesian_network.png"):
        # Generate a picture with a visualization of the network
        # path is the location where to export the picture
        graph = pydot.Dot(graph_type='digraph')

        bn_nodes = self.nodes
        for bn_node in bn_nodes:
            variable = bn_node.name
            node =  pydot.Node(variable, style="filled", fillcolor="white")
            graph.add_node(node)
            for parent in bn_node.parents:
                graph.add_edge(pydot.Edge(parent.name, node))

        graph.write_png(path)
                
                
    
    
    
        