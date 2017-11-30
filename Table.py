# A table for a CPD. Contains the combinations of variable values and parameters.
from decimal import Decimal
import copy

class Table():
    def __init__(self, variables=[], parameters=[]):
        if len(variables)>0:
            self.data = self.__buildTable(variables)
        else:
            self.data = {}
        self.__first_variable = None
    
    def __buildTable(self, variables):
        table = {}
        for combination in self.__getVariableValuesCombinations(variables):
            f_comb = frozenset(combination)
            table[f_comb] = Decimal('0.0')
        return table
    
    def __getVariableValuesCombinations(self, variables):
        # Generate all possible combinations of data
        variable_values = {}
        variable_card_prod = {}
        count = 1
        
        for variable in variables:
            variable_values[variable.name] = variable.values
            variable_card_prod[variable.name] = count
            count = count * len(variable.values)
        
        #Generate all variables combinations probabilities
        var_counters = {}
        var_counters = copy.copy(variable_card_prod)
        var_actual_values = [0] * len(variables)
        k = 0
        while (k < count):
            example = []
            for i in range(len(variables)):
                var_name = variables[i].name
                value = variable_values[var_name][var_actual_values[i]]
                var_counters[var_name] = var_counters[var_name] - 1
                if var_counters[var_name] == 0:
                    var_counters[var_name] = variable_card_prod[var_name]
                    if (var_actual_values[i] != ((len(variable_values[var_name])-1))):
                        var_actual_values[i] = var_actual_values[i] + 1
                    else:
                        var_actual_values[i] = 0
                example.append((var_name,value))
            yield example
            k+=1;
        
    def __str__(self):
        # Prints this CPD to the stdout
        result = ""    
        for example in self.data:
           str_example = ""
           feature_val_lst = []
           main_feature = ""
           for feature_val in example:
               if (self.__first_variable==None) or (feature_val[0]!=self.__first_variable.name):
                   feature_val_lst.append(feature_val[0] + "(" + str(feature_val[1]) + ")")
                   feature_val_lst.sort()
               else:
                   main_feature = feature_val[0] + "(" + str(feature_val[1]) + ")"

           if(self.__first_variable!=None):
               feature_val_lst = [main_feature] + feature_val_lst
           
           for feat_str in feature_val_lst:
               str_example += '{0:30s} '.format(feat_str)
           str_example += '{0:25s} '.format(str(self.data[example]))
           result+=str_example+"\n"
        return result
    
    # Set the parameters in the CPD
    def setParameters(self, parameters):
        # parameters is a list of parameters in the order that it appears on the table
        k=0
        for combination in self.data:
            self.data[combination]=Decimal(str(parameters[k]))
            k+=1
        
    def getParameters(self):
        return self.data
    
    def setFirstVariableToPrint(self, variable):
        self.__first_variable = variable

    def summout(self, variable, variable_values):
        # Summing out in a variable
        # Get CPD        
        table_data = self.data

        # summing out
        sum_table_data = {}
        check = {} # only to check if the work was already done
        for feature_set in table_data.keys():
            search_features = list(feature_set)
            vars = next(zip(*search_features))
            if variable in vars:
                index = vars.index(variable)
                search_features.pop(index)
                search_features = frozenset(search_features)
            if search_features not in check:
                if (len(search_features)!=len(feature_set)):
                    summ = Decimal('0.0')
                    for value in variable_values:
                        union_set = search_features.union(frozenset([(variable, value)]))
                        if union_set in table_data.keys():
                            summ += table_data[union_set]
                            check[search_features] = None

                    sum_table_data[search_features] = summ
                else:
                    sum_table_data[feature_set] = table_data[feature_set]
        
        t = Table()
        t.data = sum_table_data
        return t


    def joinFactors(self, table):       
        # Joins two factors into a single one
        CPDa = table.data
        # get features in CPDa
        CPDa_feature_list = next(zip(*list(CPDa)[0]))

        CPDb = self.data
        # get features in CPDb
        CPDb_feature_list = next(zip(*list(CPDb)[0]))

        len_a = len(CPDa_feature_list)
        len_b = len(CPDb_feature_list)
        
        if len_a == 0:
            return  CPDb

        if len_b == 0:
            return CPDa
        
        new_table_data = {}
        if (len_a !=0) or (len_b!=0):
            for feature_set_a in CPDa:
                for feature_set_b in CPDb:
                    new_feature_set = feature_set_a.union(feature_set_b)
                    check_features = next(zip(*new_feature_set))
                    if len(set(check_features)) == len(new_feature_set):
                        new_table_data[new_feature_set] = CPDa[feature_set_a] * CPDb[feature_set_b]
                   
        new_table = Table()
        new_table.data = new_table_data
        return new_table