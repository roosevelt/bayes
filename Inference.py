# Implemented inference algorithims
from decimal import Decimal
from Table import Table
import itertools as it

class Inference():
                        
    def __init__(self):
        self.any_value = "*****"                 
             
    def __buildQuery(self, query, evidences, include_any_value=True):
        # builds an understandable query to the algorithms from query and evidence
        # include_any_value determines if the returned query list and evidence list should include the values that would receive the any_value constant
        query_list = []
        ev_list = []

        if isinstance(query, str):
            if include_any_value:
                query = (query, self.any_value)
                query_list.append(query)
        else:
            if isinstance(query, tuple):
                query_list.append(query)
                
        
        if isinstance(query, list):
            for qu in query:
                if isinstance(qu, str):
                    if include_any_value:
                        query_list.append((qu, self.any_value))
                else:
                    if not include_any_value:
                        if (qu[1].strip() != self.any_value):
                            query_list.append(qu)
                    else:
                        query_list.append(qu)
                    
        for ev in evidences:
            if isinstance(ev, str):
                if include_any_value:
                    ev_list.append((ev, self.any_value))
            else:
                if (ev[1].strip() != self.any_value):
                    ev_list.append((ev))

        return query_list, ev_list    
  
    def __normalize(self, table_num, table_den):
        new_table_data ={}
        numerator = table_num.data
        denominator = table_den.data
        if isinstance(denominator, dict):
            for feature_set_d in denominator:
                for feature_set_n in numerator:
                    if feature_set_d.issubset(feature_set_n):
                        if numerator[feature_set_n] != Decimal('0.0'):
                            new_table_data[feature_set_n] = numerator[feature_set_n] / denominator[feature_set_d]
                        else:
                            new_table_data[feature_set_n]=Decimal('0.0')
        else:
            for feature_set_n in numerator:
                if numerator[feature_set_n] != Decimal('0.0'):
                    new_table_data[feature_set_n] = numerator[feature_set_n] / denominator
                else:
                    new_table_data[feature_set_n]=Decimal('0.0')
                    
        new_table = Table()
        new_table.data = new_table_data
        return new_table
    
    def __eliminationStep(self, query, bayesian_network):
        # initial CPD_set
        table_set = []
        for node in bayesian_network.nodes:
            table_set.append(node.CPD.table)
        
        # Tried to eliminate variables before, but normalization gets not so easy
        table_set = self.__eliminateComplementaryLinesFromTables(query, table_set, bayesian_network)
        node_set  = set([node.name for node in bayesian_network.nodes])
        query_vars = next(zip(*query))
        query_vars =  set(query_vars)
        query_vars = node_set.difference(query_vars)

        for variable in query_vars:
            values = bayesian_network.getNodeByName(variable).variable.values
            table_set = self.__eliminateVariableFromTables(variable, values, table_set)

        # if several tables resulted, join all of them)
        stack = table_set
        join = stack.pop(0)
        count = 0
        while len(stack) != 0:
            factor = stack.pop(0)
            
            tmp_join = join.joinFactors(factor)
            
            if len(tmp_join.data) > 0:
                join = tmp_join
                count = 0
            else:
                count += 1
                stack.append(factor)
                if (len(stack)==count):
                    break

        return join

    def __eliminateVariableFromTables(self, variable, values, table_set):
        # Eliminate a variable (as in Variable Elimination Algorithm)
        # Determine which CPDs have variable

        val_tuples = list(it.product([variable], values))

        new_table_set = []
        cpd_elim_set = []
        for table in table_set:
            features = list(table.data)[0]
            diff = features.difference(set(val_tuples))
            if len(diff) < len(features):
                cpd_elim_set.append(table)
            else:
                new_table_set.append(table)

        # Doing elimination, two steps: join and summing out factors
        join = cpd_elim_set[0]
        for i in range(len(cpd_elim_set)-1):
            join = join.joinFactors(cpd_elim_set[i+1])
            
        join = join.summout(variable, values)
        # Create new CPD set
        
        # chek if is the case of {frozenset([]): 1}, that is, all variables were eliminated
        if (len(list(join.data)[0])!=0) :
            if len(join.data)>0:
                tmp = []
                tmp.append(join)
                new_table_set.extend(tmp)

        return new_table_set
    
    def __eliminateComplementaryLinesFromTables(self, variable_list, table_set, bayesian_network):
        # Eliminate from CPD the values that don't correspond to the values of the variables
        elimination_list = []

        for variable in variable_list:                   
            for value in bayesian_network.getNodeByName(variable[0]).variable.values:
                if (value != variable[1]):
                    tupl = (variable[0],value)
                    elimination_list.append(tupl)

        elim_set = frozenset(elimination_list)

        new_table_array = []        
        for table in table_set:
            new_table_data = {}
            for feature in table.data:
                if not(len(feature.difference(elim_set)) < len(feature)):
                    new_table_data[feature] = table.data[feature]
            t = Table()
            t.data = new_table_data
            new_table_array.append(t)

        return new_table_array
    
    
    def variableElimination(self, query, evidence, bayesian_network):
        # Returns exact inference made using Variable Elimination Algorithm

        # Adapt to different types of entries
        query_list, ev_list = self.__buildQuery(query, evidence, False)
        variables_with_values = query_list + ev_list

        # Calculate
        if len(query_list) ==1:
            query_vars = [query_list[0][0]]
            query_values = [query_list[0][1]]
        else:
            lst = list(zip(*query_list))
            query_vars =  list(lst[0])
            query_values = list(lst[1])

        ev_vars = []
        ev_values = []
        if len(ev_list) != 0:
            if len(ev_list) ==1:
                ev_vars = [ev_list[0][0]]
                ev_values = [ev_list[0][1]]
            else:
                lst = list(zip(*ev_list))
                ev_vars =  list(lst[0])
                ev_values = list(lst[1])

        set_query_ev_vars = set(query_vars+ev_vars)
        query_ev_vars = query_vars + ev_vars
        query_ev_values = query_values + ev_values
        if len(set_query_ev_vars) != len(query_ev_vars):
            for k in range(len(query_ev_vars)):
                for w in range(k, len(query_ev_vars)):
                    if query_ev_vars[k]==query_ev_vars[w]:
                        if query_ev_values[k] != query_ev_values[w]:
                            return Table()
        
        
        num = self.__eliminationStep(variables_with_values, bayesian_network)
        # Normalize
        table = num
        if len(ev_list) != 0:
            den = self.__eliminationStep(ev_list, bayesian_network)
            table = self.__normalize(num, den)
        tmp = []
        tmp.append(table)
        table_set = self.__eliminateComplementaryLinesFromTables(query_list, tmp, bayesian_network)

        return table_set[0]