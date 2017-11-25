# Oriented Edge of a Bayesian Network

class Edge():
    def __init__(self, node_origin, node_destiny):
        self.origin = node_origin
        self.destiny = node_destiny
        
    def __str__(self):
        buffer = "Edge from node " + self.origin.name + " to node " + self.destiny.name
        return buffer
        
        