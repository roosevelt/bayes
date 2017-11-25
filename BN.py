# Bayesian network class 
from Node import Node
from Edge import Edge

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
            
                
                
    
    
    
        