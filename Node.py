# A Bayesian Node, with its properties and associated CPD
class Node():
    
    def __init__(self, name, CPD):
        self.name = name
        self.CPD = CPD
        self.parents = []
        
    
        