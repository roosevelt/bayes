# Bayesian network class 
from Variable import Variable
from CPD import CPD
from Node import Node
from Edge import Edge
import pydot
from decimal import Decimal
import copy

class BayesianNetwork:
    
    def __init__(self, nodes=[]):
        self.nodes = nodes
        if (len(self.nodes)>0):
            self.edges = self.__buildEdges(self.nodes)
        else:
            self.edges = []
            
        self.__node_dict = self.__buildNodeDict(self.nodes)
        
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
        
    def loadBIF(self, input_path):
        content = open(input_path, "r", encoding='utf-8').read();
        
        lines = content.splitlines()
                        
        # First part: parse
        variables = []
        parents = {}
        possible_values = {}
        cpd = {}
        getCPD = False
        for line in lines:
            # check if variable, if it is structure
            if line.strip().startswith("variable"):
                variable_name = line.split(" ")[1]
                variable_name.replace("{", "").strip()
                if variable_name not in variables:
                    variables.append(variable_name)
                    
            # check if parents
            if line.strip().startswith("probability"):
                variable_parents_children = line.split("(")[1].split(")")[0].strip()
                if "|" in variable_parents_children:
                    variable_parents =  variable_parents_children.split("|")[1].strip().replace(" ", "").split(",")
                    variable = variable_parents_children.split("|")[0].strip()
                else:
                    variable_parents = []
                    variable = variable_parents_children.strip()
                    
                parents[variable] = variable_parents
                cpd_variable = variable
                getCPD =True
                
            # get values of variables
            if line.strip().startswith("type"):
                values = line.strip().split("{")[1].strip().split("}")[0].replace(" ", "").split(",")
                possible_values[variable_name]=values
                 
            # now cpd
            if getCPD == True and not(line.strip().startswith("probability")):
                # When the brackets are closed, end parsing a cpd
                if "}" in line:
                    getCPD = False
                else:
                    if line.strip().startswith("table"):
                        probs = line.strip().replace("table", "").replace(" ", "").replace(";", "").split(",")
                        k=0
                        for value in possible_values[cpd_variable]:
                            tupl_cpd_var = (cpd_variable, value)
                            if cpd_variable not in cpd:
                                cpd[cpd_variable] = {}
                                cpd[cpd_variable][frozenset(set([tupl_cpd_var]))] = Decimal(probs[k]) 
                            else:
                                cpd[cpd_variable][frozenset(set([tupl_cpd_var]))] = Decimal(probs[k])
                            k+=1
                            
                        # TODO
                    else:
                        values = line.strip().split("(")[1].split(")")[0].strip()
                        if "," in values:                   
                            values = values.replace(" ","").strip().split(",")
                        else:
                            values = [values]
                
                        probs = line.strip().split(")")[1].replace(" ", "").replace(";", "").split(",")
                    
                        k=0
                        arr_keys = []
                        for parent in parents[cpd_variable]:
                            tupl = (parent, values[k])
                            arr_keys.append(tupl)
                            k+=1
                            
                        k = 0    
                        for value in possible_values[cpd_variable]:
                            arr = copy.deepcopy(arr_keys)
                            tupl_cpd_var = (cpd_variable, value)
                            arr.insert(0, tupl_cpd_var)
                            if cpd_variable not in cpd:
                                cpd[cpd_variable] = {}
                                cpd[cpd_variable][frozenset(set(arr))] = Decimal(probs[k]) 
                            else:
                                cpd[cpd_variable][frozenset(set(arr))] = Decimal(probs[k]) 
                            k+=1
          
        # Second part: build bn
        var_dict = {}
        for var_name in variables:
            variable = Variable(var_name, possible_values[var_name])
            var_dict[var_name]=variable
            
        node_list = []
        for var_name in variables:
            variable = var_dict[var_name]
            condition_variables = []
            for parent in parents[var_name]:
                condition_variables.append(var_dict[parent])                                
            cpd_obj = CPD(variable, condition_variables)
            cpd_obj.table = cpd[var_name]
            node = Node(cpd_obj)
            node_list.append(node)

        self.nodes = node_list
        self.edges = self.__buildEdges(self.nodes)
        self.__node_dict = self.__buildNodeDict(self.nodes)
        
        
            
        
    def getNodeByName(self, name):
        return self.__node_dict[name]
        
    
    def __buildNodeDict(self, nodes):
        node_dict = {}
        for node in nodes:
            node_dict[node.name] = node
        return node_dict

    
                
                
    
    
    
        