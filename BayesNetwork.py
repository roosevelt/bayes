# -*- coding: iso-8859-1 -*-
import string
import array
import os.path
from storage import Storage
import re
import random
import numpy as np
import copy
import math
import pydot
import time
import timer
import itertools as it
import sys
import os
import decimal
import pp
import imp
import networkx
import collections
import gc
import operator
from scipy.stats import chi2
from scipy.stats import ttest_rel
sys.path.append("./lib/progressbar/")
from progressbar import ProgressBar


# author: Roosevelt de Lima Sardinha

ppservers = ()
job_server = pp.Server(ppservers=ppservers, socket_timeout=36000)

class BayesNetwork():
    def __init__(self):
        self.BayesianNetwork = Storage()
        self.any_value = "*****"
        self.showProgress = True
        self.showTime = True
        self.cores = 2
        self.smoothing_param = decimal.Decimal('0.0000001')
        self.very_small_value = decimal.Decimal('0.0000001')
        self.adjacencies_dict = None
        decimal.getcontext().prec = 10
        
        
    def parse(self, input_path):
        # parses a bayesian network
        extension =  input_path.strip().split(".")[-1].strip()
        if extension == "eg":
            bn = self.parse_encog(input_path)
        elif extension == "bif":
            bn = self.parse_bif(input_path)
        return bn
        
        
    
    def parse_bif(self, input_path):
        content = unicode(open(input_path, "r").read(), "iso-8859-1");
        
        lines = content.splitlines()
        
        BayesianNetwork = Storage()
                
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
                                cpd[cpd_variable][frozenset(set([tupl_cpd_var]))] = decimal.Decimal(probs[k]) 
                            else:
                                cpd[cpd_variable][frozenset(set([tupl_cpd_var]))] = decimal.Decimal(probs[k])
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
                                cpd[cpd_variable][frozenset(set(arr))] = decimal.Decimal(probs[k]) 
                            else:
                                cpd[cpd_variable][frozenset(set(arr))] = decimal.Decimal(probs[k]) 
                            k+=1
          
        # Second part: build bn
        for node in variables:
            BayesianNetwork[node] = Storage()
            BayesianNetwork[node].CPD = cpd[node]
            BayesianNetwork[node].Values = possible_values[node]
            BayesianNetwork[node].Parents = parents[node]
            
        self.BayesianNetwork = BayesianNetwork
        return BayesianNetwork
                
        
    def parse_encog(self,input_path):
        #input_path = "../data/EncogBayesData/alarm.eg"

        content = unicode(open(input_path, "r").read(), "iso-8859-1");

        lines = content.splitlines()

        BayesianNetwork = Storage()

        # Each line to parse
        for line in lines:
            # Line of the structure
            if (line.strip().startswith("contents=")):
                probs = re.findall("P\(.+?\)", line)
                for prob in probs:
                    div = prob.split("|")
                    if (len(div) >= 1):
                        # Get the name of the node
                        node = div[0].split("[")[0].split("(")[1]
                        BayesianNetwork[node] = Storage()

                        BayesianNetwork[node].Parents = []
                        if len(div)==2:
                            evidences = div[1].replace(")", "").strip()
                            # Get evidences
                            evs = evidences.split(",")
                            for ev in evs:
                                BayesianNetwork[node].Parents.append(ev.strip())

                        #Get values
                        vals = div[0].split("[")
                        vals = vals[1].split("]")[0]
                        vals = vals.split(",")
                        BayesianNetwork[node].Values = []
                        for val in vals:
                            BayesianNetwork[node].Values.append(val)

            # CPDs
            if (line.strip().startswith("P(")):
                line = line.replace("P(","")
                line = line.split("|")
                node_val = line[0].strip()
                node = line[0].split("=")[0].replace("+","").replace("-","").replace(")", "").strip();

                features_set = []
                if BayesianNetwork[node].CPD is None:
                    BayesianNetwork[node].CPD = {}
                if len(line) == 2:
                    evidences = line[1].split(")")[0].split(",")
                    for ev in evidences:
                        if ev.strip().startswith("+"):
                            if len(BayesianNetwork[ev.replace("+","").strip()].Values) == 2:
                                features_set.append((ev.replace("+","").strip(), BayesianNetwork[ev.replace("+","").strip()].Values[0]))
                            else:
                                features_set.append((ev.replace("+","").strip(), True))

                        elif ev.strip().startswith("-"):
                            if len(BayesianNetwork[ev.replace("-","").strip()].Values) == 2:
                                features_set.append((ev.replace("-","").strip(), BayesianNetwork[ev.replace("-","").strip()].Values[1]))
                            else:
                                features_set.append((ev.replace("-","").strip(), False))
                        elif ev.find("=") != -1:
                            ev = ev.split("=")
                            feature = ev[0]
                            value = ev[1]
                            features_set.append((feature.strip(), value))
                        else:
                            print "--------------Invalid Feature!!-------------"

                node_val = node_val.split(")")[0]
                if node_val.find("=") != -1:
                    feat_and_val = node_val.split("=")
                    feature = feat_and_val[0]
                    if not feat_and_val[1].isdigit():
                        value = feat_and_val[1]
                    else:
                        value = BayesianNetwork[feature].Values[int(feat_and_val[1])]

                if node_val.strip().startswith("+"):
                    feature = node_val.strip().replace("+", "")
                    if len(BayesianNetwork[feature].Values) == 2:
                        value = BayesianNetwork[feature].Values[0]
                    else:
                        value = True

                if node_val.strip().startswith("-"):
                    feature = node_val.strip().replace("-", "")
                    if len(BayesianNetwork[feature].Values) == 2:
                        value = BayesianNetwork[feature].Values[1]
                    else:
                        value = False


                features_set.append((feature.strip(), value))

                # Get the probability of the features to happen
                if len(line)==2:
                    prob_str = line[1].split("=")[-1].strip()
                elif len(line)==1:
                    prob_str = line[0].split("=")[-1].strip()

                if prob_str=='None':
                    probability = None
                else:
                    probability = decimal.Decimal(prob_str)
                BayesianNetwork[node].CPD[frozenset(features_set)] = probability

        self.BayesianNetwork = BayesianNetwork
        return BayesianNetwork

    def export_to_bnt(self, function_name,  output, struct_only=False, bn=None):
        # Export bayes network to the format used by BNT
        buffer = ""
        
        if bn is None:
            bn = self.BayesianNetwork
        
        # Get variables in topological order; needed by BNT
        variables = self.getTopologicalOrder(bn)
        # Get number of variables
        N = len(variables)
        
        # Start building function
        buffer+="function bnet = "+ function_name +"()\n"
        buffer+="N="+ str(N) + ";\n"
        buffer+="dag = zeros("+str(N)+","+str(N)+");\n"
        
        # Build graph structure
        for variable in variables:
            for parent in bn[variable].Parents:
                buffer+="dag("+ str(variables.index(parent)+1)  + ","+ str(variables.index(variable)+1) +")=1;\n"
        
        # Build the data about the size of the nodes
        buffer+="node_sizes= zeros(1,"+ str(N) +");\n"
        for variable in variables:
            buffer+="node_sizes("+str(variables.index(variable)+1)+")="+str(len(bn[variable].Values))+";\n"
        
        # Build bayesian network
        buffer+= "bnet = mk_bnet(dag, node_sizes);\n"
        
        # Define the CPDs if struct_only option is disabled
        if not struct_only:
            for variable in variables:
                parents = sorted(bn[variable].Parents, key=lambda str: variables.index(str))
                buffer+="CPT= zeros("
                for parent in parents:
                    buffer+= str(len(bn[parent].Values)) + ", "
                    
                if len(parents)!=0:
                    buffer+= str(len(bn[variable].Values)) + ");\n" 
                else:
                    buffer+= str(len(bn[variable].Values)) + ",1);\n"
                
                for line in bn[variable].CPD:
                    dict ={}
                    for tpl in line:
                        dict[tpl[0]] = tpl[1];
                    buffer+="CPT("
                    for parent in parents:
                        val =  dict[parent]
                        buffer+= str(bn[parent].Values.index(val) + 1) + ","
                    val =  dict[variable]
                    buffer+=str(bn[variable].Values.index(val) + 1) + ")=" + str(bn[variable].CPD[line]) +";\n"
                
                buffer+="CPT=reshape(CPT, 1, prod(size(CPT)));\n"
                buffer+="bnet.CPD{"+str(variables.index(variable)+1)+"} = tabular_CPD(bnet, "+str(variables.index(variable)+1)+", 'CPT', CPT);\n"
                #buffer+="bnet.CPD{"+str(variables.index(variable)+1)+"}=CPT;\n"

        out = open(output, "w")
        out.write(buffer)
        out.close()

    def export_encog(self, output, bn=None):
        # Exports a bayesian network to the encog file format, for result persistence
        if bn is None:
            bn=self.BayesianNetwork
        out = open(output, "w")

        buffer = ""
        buffer +="[BAYES-NETWORK]\n"
        buffer +="[BAYES-NETWORK:BAYES-PARAM]\n"
        buffer +="queryType=EnumerationQuery\n"
        buffer +="query=P(|)\n"
        buffer +="contents="

        # Build contents line
        contents = ""
        for node in bn:
            if len(bn[node].Parents)!=0:
                contents += "P(" + node + "[" + ",".join(bn[node].Values)  + "]|" +  ",".join(bn[node].Parents) + ") "
            else:
                contents += "P(" + node + "[" + ",".join(bn[node].Values)  + "]) "
                
        buffer+=contents + "\n"
        buffer+="[BAYES-NETWORK:BAYES-PROPERTIES]\n"
        buffer+="[BAYES-NETWORK:BAYES-TABLE]\n"


        for node in bn:
            cpd_contents = ""
            for comb in bn[node].CPD:
                ev = []
                for tpl in comb:
                    if tpl[0]==node:
                        node_text = tpl[0]+ "=" + str(tpl[1])
                    else:
                        ev.append(tpl[0]+ "=" + str(tpl[1]))
                
                if len(ev)>0:
                    cpd_contents += "P(" + node_text + "|" + ",".join(ev) + ")=" + str(bn[node].CPD[comb]) + "\n"
                else:
                    cpd_contents += "P(" + node_text + ")=" + str(bn[node].CPD[comb]) + "\n"
   
            buffer +=cpd_contents

       
        out.write(buffer)
        out.close()


    def printCPD(self, node):
        # Prints the CPD of a node in the terminal for evaluation
        if isinstance(node, str):
            cpd = self.BayesianNetwork[node].CPD
        else:
            cpd = node

        for example in cpd:
           str_example = ""
           feature_val_lst = []

           for feature_val in example:
               feature_val_lst.append(feature_val[0] + "(" + str(feature_val[1]) + ")")
               feature_val_lst.sort()

           for feat_str in feature_val_lst:
               str_example += '{0:30s} '.format(feat_str)
           str_example += '{0:25s} '.format(str(cpd[example]))

           print str_example

    def getProb(self, query, evidence=[], cpd=None):
        # Returns the probability from a CPD
        if cpd is None:
            cpd = self.BayesianNetwork[query[0]].CPD
        tmp = []
        tmp.append(query)
        ev_query = tmp + evidence
        key = frozenset(ev_query)
        if cpd is None:
            prob = self.BayesianNetwork[query[0]].CPD[key]
        else:
            prob = cpd[key]
        return prob

    def joinFactors(self, factorA, factorB, bn=None):       
        if bn is None:
            bn=self.BayesianNetwork

        # Joins two factors into a single one
        if isinstance(factorA, str):
            CPDa = bn[factorA].CPD
        else:
            CPDa = factorA

        # get features in CPDa
        CPDa_feature_list = zip(*CPDa.keys()[0])[0]

        if isinstance(factorB, str):
            CPDb = bn[factorB].CPD
        else:
            CPDb = factorB

        # get features in CPDb
        CPDb_feature_list = zip(*CPDb.keys()[0])[0]

        newCPD = {}

        len_a = len(CPDa_feature_list)
        len_b = len(CPDb_feature_list)
        if (len_a !=0) or (len_b!=0):
            for feature_set_a in CPDa:
                for feature_set_b in CPDb:
                    new_feature_set = feature_set_a.union(feature_set_b)
                    check_features = zip(*new_feature_set)[0]
                    if len(set(check_features)) == len(new_feature_set):
                        newCPD[new_feature_set] = CPDa[feature_set_a] * CPDb[feature_set_b]

        if len_a == 0:
            new_CPD = CPDb

        if len_b == 0:
            new_CPD = CPDa

        return newCPD

    def summout(self, node, variable, bn=None):
        # Summing out in a variable
        # Get CPD
        if bn is None:
            bn = self.BayesianNetwork
        
        if isinstance(node, str):
            CPD = bn[node]['CPD']
        else:
            CPD = node

        # summing out
        newCPD = {}
        check = {} # only to check if the work was already done
        for feature_set in CPD.keys():
            search_features = list(feature_set)
            vars = zip(*search_features)[0]
            if variable in vars:
                index = vars.index(variable)
                search_features.pop(index)
                search_features = frozenset(search_features)
            if search_features not in check:
                if (len(search_features)!=len(feature_set)):
                    summ = decimal.Decimal('0.0')
                    values = bn[variable]['Values']
                    for value in values:
                        union_set = search_features.union(frozenset([(variable, value)]))
                        if union_set in CPD.keys():
                            summ += CPD[union_set]
                            check[search_features] = None

                    newCPD[search_features] = summ
                else:
                    newCPD[feature_set] = CPD[feature_set]

        return newCPD

    def eliminateVariable(self, variable, CPD_set, bn=None):
        # Eliminate a variable (as in Variable Elimination Algorithm)
        # Determine which CPDs have variable
        if bn is None:
            bn = self.BayesianNetwork

        val_tuples = list(it.product([variable], bn[variable]['Values']))

        newCPD_set = []
        cpd_elim_set = []
        for cpd in CPD_set:
            features = cpd.keys()[0]
            diff = features.difference(set(val_tuples))
            if len(diff) < len(features):
                cpd_elim_set.append(cpd)
            else:
                newCPD_set.append(cpd)

        # Doing elimination, two steps: join and summing out factors
        join = cpd_elim_set[0]
        for i in range(len(cpd_elim_set)-1):
            join = self.joinFactors(join, cpd_elim_set[i+1], bn)
            
        join = self.summout(join, variable, bn)
        # Create new CPD set
        
        # chek if is the case of {frozenset([]): 1}, that is, all variables were eliminated
        if (len(join.keys()[0])!=0) :
            if len(join)>0:
                tmp = []
                tmp.append(join)
                newCPD_set.extend(tmp)

        return newCPD_set

    def normalize(self, numerator, denominator):
        newCPD ={}
        if isinstance(denominator, dict):
            for feature_set_d in denominator:
                for feature_set_n in numerator:
                    if feature_set_d.issubset(feature_set_n):
                        if numerator[feature_set_n] != decimal.Decimal('0.0'):
                            newCPD[feature_set_n] = numerator[feature_set_n] / denominator[feature_set_d]
                        else:
                            newCPD[feature_set_n]=decimal.Decimal('0.0')
        else:
            for feature_set_n in numerator:
                if numerator[feature_set_n] != decimal.Decimal('0.0'):
                    newCPD[feature_set_n] = numerator[feature_set_n] / denominator
                else:
                    newCPD[feature_set_n]=decimal.Decimal('0.0')
        return newCPD

    def buildQuery(self, query, evidences, include_any_value=True):
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

    def eliminateComplementaryLines(self, variable_list, CPD_set, bn=None):
        # Eliminate from CPD the values that don't correspond to the values of the variables
        if bn is None:
            bn = self.BayesianNetwork       
        
        elimination_list = []

        for variable in variable_list:                   
            for value in bn[variable[0]]['Values']:
                if (value != variable[1]):
                    tupl = (variable[0],value)
                    elimination_list.append(tupl)

        elim_set = frozenset(elimination_list)

        newCPD_set = []        
        for cpd in CPD_set:
            newCPD = {}
            for feature in cpd:
                if not(len(feature.difference(elim_set)) < len(feature)):
                    newCPD[feature] = cpd[feature]
            newCPD_set.append(newCPD)

        return newCPD_set


    def __eliminationStep(self, query, bn=None):
        if bn is None:
            bn = self.BayesianNetwork

        # initial CPD_set
        CPD_set = []
        for node in bn:
            CPD_set.append(bn[node]['CPD'])
        
        # Tried to eliminate variables before, but normalization gets not so easy
        CPD_set = self.eliminateComplementaryLines(query,CPD_set, bn)
        node_set  = set(bn.keys())
        query_vars = zip(*query)[0]
        query_vars =  set(query_vars)
        query_vars = node_set.difference(query_vars)

        for variable in query_vars:
            CPD_set = self.eliminateVariable(variable, CPD_set, bn)

        # if several tables resulted, join all of them

        stack = CPD_set
        join = stack.pop(0)
        count = 0
        while len(stack) != 0:
            factor = stack.pop(0)
            
            tmp_join = self.joinFactors(join, factor, bn)
            
            if len(tmp_join) > 0:
                join = tmp_join
                count = 0
            else:
                count += 1
                stack.append(factor)
                if (len(stack)==count):
                    break

        return join;

    def variableEliminationInference(self, query, evidence=[], bn=None):
        # Returns exact inference made using Variable Elimination Algorithm

        if bn is None:
            bn = self.BayesianNetwork

        # Adapt to different types of entries
        query_list, ev_list = self.buildQuery(query,evidence, False)
        variables_with_values = query_list + ev_list

        # Calculate
        if len(query_list) ==1:
            query_vars = [query_list[0][0]]
            query_values = [query_list[0][1]]
        else:
            lst = zip(*query_list)
            query_vars =  list(lst[0])
            query_values = list(lst[1])

        ev_vars = []
        ev_values = []
        if len(ev_list) != 0:
            if len(ev_list) ==1:
                ev_vars = [ev_list[0][0]]
                ev_values = [ev_list[0][1]]
            else:
                lst = zip(*ev_list)
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
                            return {}
        
        num = self.__eliminationStep(variables_with_values, bn)
        # Normalize
        CPD = num
        if len(ev_list) != 0:
            den = self.__eliminationStep(ev_list, bn)
            CPD = self.normalize(num, den)
        tmp = []
        tmp.append(CPD)
        CPD_set = self.eliminateComplementaryLines(query_list,tmp, bn)

        return CPD_set[0]


    def getMarkovBlanket(self, node, bn=None):
        # Returns a list with the name of the nodes in the markov blaket of a node
        if bn is None:
            bn = self.BayesianNetwork
        mb = set(bn[node].Parents)

        for variable in bn:
            if variable != node:
                if node in bn[variable].Parents:
                    mb = mb.union(bn[variable].Parents)
                    mb = mb.union(set([variable]))
        mb = list(mb)

        if node in mb:
            mb.remove(node)
        return mb

    def simpleSampling(self, values, probs):
        # generates a random number and selects a value according to its probability
        rnd = decimal.Decimal(random.random())
        summ = 0
        for i in range(len(values)):
            if (rnd >= summ) and (rnd < (summ + probs[i])):
                result = values[i]
            summ += probs[i]
        return result

    def gibbsSampler(self, variables, initial_state ,iterations):
        # Gibbs Sampler
        # Parameters :
        # - The variables to sample
        # - The initial state of the variables to sample
        # - No of iteractions
        bn  =  self.BayesianNetwork

        samples_matrix = np.array([initial_state],dtype="|S50")
        for k in range(iterations):
            for i in range(len(variables)):
                # Markov blanket of the node
                mb = self.getMarkovBlanket(variables[i])

                # Build list of pairs: evidence -> value
                evidences = list(set(variables).intersection(set(mb)))
                if variables[i] in evidences:
                    evidences.remove(variables[i])

                ev_tupls = []
                for ev in evidences:
                    index = variables.index(ev)
                    a = samples_matrix[(samples_matrix.shape[0]-1), index]
                    ev_tupls.append((ev, a))

                # Calculates conditional probability using exact inference
                cpd = self.variableEliminationInference(variables[i], ev_tupls)
                var_values = bn[variables[i]].Values
                probabilities = []
                for value in var_values:
                    probabilities.append(self.getProb((variables[i], value), ev_tupls, cpd))

                # New value for the variable
                value = self.simpleSampling(var_values, probabilities)

                # Update sample
                new_sample = np.array([])
                new_sample = samples_matrix[samples_matrix.shape[0]-1]
                new_sample[i] = value

                # Stores in sample history
                samples_matrix = np.vstack((samples_matrix, new_sample))

        return samples_matrix
    
    def gibbsSamplingInference(self, query, evidence, bn=None):
        # Answers an inference query. Ex.: P(B | C=True, M=True)
        if bn is None:
            bn = self.BayesianNetwork
            
        query_list, ev_list = self.buildQuery(query,evidence)
        alg_query = query_list + ev_list
        iterations = 1000
        # Get variables and initial_state for gibbs sampling
        variables = []
        initial_state = []
        for tupl in alg_query:
            variables.append(tupl[0])
            initial_state.append(bn[tupl[0]].Values[0])

        sample_matrix = self.gibbsSampler(variables, initial_state, iterations)

        filtered_matrix = sample_matrix
        for i in range(len(alg_query)):
            if (i != 0):
                ev = alg_query[i][0]
                value = alg_query[i][1]
                var_position = variables.index(ev)
                filtered_matrix = filtered_matrix[(filtered_matrix[:, var_position]==value)]

        total = filtered_matrix.shape[0]
        query_var = alg_query[0][0]
        var_position = variables.index(query_var)
        cpd = {}
        for value in bn[query_var].Values:
            idx = (filtered_matrix[:, var_position]==str(value))
            count_matrix = filtered_matrix[idx]
            count = count_matrix.shape[0]
            probability = decimal.Decimal(count)/decimal.Decimal(total)
            example = list(alg_query)
            example[0] = (query_var, value)
            cpd[frozenset(example)] = probability

        tmp = []
        tmp.append(cpd)
        tmp2 = []
        tmp2.append(alg_query[0])

        cpd_set = self.eliminateComplementaryLines(tmp2, tmp, bn)
        return cpd_set[0]

    def generateVariableValuesCombinations(self, variables, bn=None):
        # Generate all possible combinations of data
        if bn is None:
            bn = self.BayesianNetwork
        variable_values = {}
        variable_card_prod = {}
        count = 1
        variables.reverse()
        for variable in variables:
            variable_values[variable] = bn[variable]['Values']
            variable_card_prod[variable] = count
            count = count * len(bn[variable]['Values'])
        variables.reverse()

        #Generate all variables combinations probabilities
        var_counters = {}
        var_counters = copy.copy(variable_card_prod)
        var_actual_values = [0] * len(variables)
        examples = []
        k = 0
        while (k < count):
            example = []
            for i in range(len(variables)):
                var_name = variables[i]
                value = variable_values[var_name][var_actual_values[i]]
                var_counters[var_name] = var_counters[var_name] - 1
                if var_counters[var_name] == 0:
                    var_counters[var_name] = variable_card_prod[var_name]
                    if (var_actual_values[i] != ((len(variable_values[var_name])-1))):
                        var_actual_values[i] = var_actual_values[i] + 1
                    else:
                        var_actual_values[i] = 0
                example.append((var_name.encode('iso-8859-1'),value.encode('iso-8859-1')))
            yield example
            k+=1;


    def generateDataProbabilityTable(self, bn=None):
        # Generate the probabilities of all possible examples
        if bn is None:
            bn = self.BayesianNetwork
        variables = bn.keys()

        join = bn[variables[0]].CPD
        for i in range(1, len(variables)):
            join = self.joinFactors(join, bn[variables[i]].CPD)

        return join


    def generateDataProbTableFromDataFreq(self, variables,data, bn=None):
        # Generate the probabilities of all possible examples
        if bn is None:
            bn = self.BayesianNetwork
        prob_table = {}
        for comb in self.generateVariableValuesCombinations(variables, bn):
            prob_table[frozenset(comb)] = decimal.Decimal('0.0');
        n = decimal.Decimal(data.shape[0])
        # Compute the probabilities
        for example in data:
            tmp = []
            for k in range(len(variables)):
                tmp.append((variables[k], example[k]))
            prob_table[frozenset(tmp)] = prob_table[frozenset(tmp)] + (decimal.Decimal('1.0')/n)
        # Return a CPD
        return prob_table


    def getTopologicalOrder(self, bn=None):
        if bn is None:
            bn = self.BayesianNetwork

        variables = bn.keys()

        # Generate children list
        children = {}
        for node in variables:
            children[node] = []

        for node in variables:
            for parent in bn[node].Parents:
                children[parent].append(node)

        # Generates a list in a topological order for the network
        to_visit = []
        for node in variables:
            if len(bn[node].Parents) == 0:
                to_visit.append(node)

        node_order = []
        while len(to_visit)!=0:
            actual_node = to_visit.pop(0)
            node_order.append(actual_node)
            for child in children[actual_node]:
                tmp = node_order+to_visit
                if child not in tmp:
                    to_visit.append(child)

        return node_order
        
    def verbose2bnt(self, input, output):
        # Converts data in the verbose format to the bnt format.
        bn=self.BayesianNetwork
        
        variables, data = self.loadDataset(input)
        
        out = open(output, 'w')
        out.write(",".join(variables) + "\n") 
       
        for example in data:
            example_bnt = []
            for k in range(len(variables)):
                val = example[k]
                var = variables[k]
                if val.strip()!="":
                    example_bnt.append(str(bn[var].Values.index(val) + 1))
                else:
                    example_bnt.append("")
            out.write(",".join(example_bnt) + "\n")    
        out.close()
        
    def forwardSampling(self, num_samples, out_file, missing_prob={}, type='verbose'):
        # Generates a dataset for a large BN by saving time and memory
        # Type can be verbose (write full name of the value) or bnt (matlab bnt compatible format)
        bn = self.BayesianNetwork
        variables = bn.keys()

        # Missing data entry
        if isinstance(missing_prob, int) or isinstance(missing_prob, float) or isinstance(missing_prob, decimal.Decimal):
            num = missing_prob
            missing_prob = {}
            for variable in variables:
                missing_prob[variable] = num

        # Order of choosing nodes (Topological order)
        node_order = self.getTopologicalOrder(bn)

        # Open file to receive results
        out = open(out_file, "w")
        out.write(",".join(node_order) + "\n")
        # Now the sampling
        k = 0
        while k < num_samples:
            example = []
            result = None
            tpls = []
            for node in node_order:
                CPD = bn[node].CPD
                if result is not None:
                    CPD = self.eliminateComplementaryLines(tpls, [CPD], bn)[0]
                result = self.simpleSampling(CPD.keys(), CPD.values())
                tmp = list(result)
                if  len(result) ==1:
                    tmp = tmp[0]
                    var = tmp[0]
                    value = tmp[1]
                else:
                    tmp = zip(*tmp)
                    vars = tmp[0]
                    values = tmp[1]
                    index = vars.index(node)
                    var = vars[index]
                    value = values[index]

                if variable in missing_prob:
                    miss = missing_prob[variable]
                    miss_result = self.simpleSampling([True, False], [miss, 1-miss])
                else:
                    miss_result = False

                if not miss_result:
                    if type=='verbose':
                        example.append(value)
                    elif type=='bnt':
                        value_index = bn[var].Values.index(value) + 1
                        example.append(str(value_index))                        
                else:
                    example.append("")
                tpls.append((var, value))

            out.write(",".join(example) + "\n")
            k+=1
        out.close()


    def generateDatasetFromBN(self, num_samples, out_file, missing_prob={}):
        ## Generates an artificial dataset with the distribution in the CPDs
        ## Generates complete datasets or datasets with missing data
        bn = self.BayesianNetwork
        variables = bn.keys()

        # Missing data entry
        if isinstance(missing_prob, int) or isinstance(missing_prob, float) or isinstance(missing_prob, decimal.Decimal):
            for variable in variables:
                missing_prob[variable] = missing_prob

        prob_table = self.generateDataProbabilityTable()

        #Generate examples by sampling from the combinations
        keys = prob_table.keys()
        values = []
        for key in keys:
            values.append(prob_table[key])

        # Formatting
        to_file = []
        to_file.append(",".join(variables))
        for i in range(num_samples):
            tmp =[]
            sampling = self.simpleSampling(keys, values)
            dic = {}
            for el in sampling:
                var_name = el[0]
                var_value = el[1]
                dic[var_name] = var_value

            for variable in variables:
                if variable in missing_prob:
                    miss = missing_prob[variable]
                    miss_result = self.simpleSampling([True, False], [miss, 1-miss])
                else:
                    miss_result = False

                if not miss_result:
                    tmp.append(dic[variable])
                else:
                    tmp.append("")
                miss_result = False

            to_file.append(",".join(tmp))
        out = open(out_file, "w")
        out.write("\n".join(to_file))
        out.close()

    def loadDataset(self, dataset):
        print "Loading Data..."
        input_file = open(dataset, "r")
        content = input_file.read()
        content_list = content.split("\n")
        for k in range(len(content_list)):
            line = content_list[k]
            if line.strip()!="":
                line_list = line.split(",")
                new_sample = np.array(line_list, dtype="|S50")
                if k!=0:
                    if k!=1:
                        samples = np.vstack((samples, new_sample))
                    else:
                        samples = new_sample
                else:
                    variables = list(new_sample)
                    
        print "Number of variables: " + str(len(variables))
        print "Number of examples: " + str(samples.shape[0])
        miss = self.countMissing(variables, samples)
        print str(miss * 100) + "% of examples missing some data."
        print "Finished loading data."
        print ""
        return variables, samples

    def randomParameterBNInitialization(self, bn=None, nodes_to_consider=[], overwrite=False):
        # Set random parameters to a network that has no parameters    
        if (len(nodes_to_consider) == 0):
            nodes_to_consider = bn.keys()
            
        for node in nodes_to_consider:
            # In case we don't have CPDs but need to initialize
            if ("CPD" not in bn[node]) or overwrite:
                bn[node].CPD ={}
                for comb in self.generateVariableValuesCombinations([node]+bn[node].Parents, bn):
                    f_comb = frozenset(comb)
                    bn[node].CPD[f_comb] = decimal.Decimal('0.0')

            CPD = bn[node].CPD
            parents = bn[node].Parents
            
            node_possibilities=[]
            for node_possibility in self.generateVariableValuesCombinations([node], bn):
                node_possibilities.append(node_possibility) 
                
            for parents_values in self.generateVariableValuesCombinations(parents, bn):
                counter = 0
                total = 1
                for node_value in node_possibilities:
                    comb  = node_value + parents_values
                    f_comb = frozenset(comb)
                    random_num = decimal.Decimal(random.random())
                    if counter == 0:
                        CPD[f_comb] = random_num
                        total  = total -  random_num
                    elif counter == (len(node_possibilities) -1):
                        CPD[f_comb] = total
                    else:
                        CPD[f_comb] = random_num * total
                        total =  total - (random_num * total)
                    counter+=1
        return bn

    def fixedParameterBNInitialization(self, parameter_value, bn=None, nodes_to_consider=[], overwrite=False):
        # Set fixed parameters to a network that has no parameters
        if bn is None:
            bn = self.BayesianNetwork

        if len(nodes_to_consider) ==0:
            nodes_to_consider = bn.keys()
            
        for node in nodes_to_consider:
            # In case we don't have CPDs but need to initialize
            if ("CPD" not in bn[node]) or (overwrite):
                bn[node].CPD ={}
                for comb in self.generateVariableValuesCombinations([node]+bn[node].Parents, bn):
                    f_comb = frozenset(comb)
                    bn[node].CPD[f_comb] = decimal.Decimal('0.0')

            CPD = bn[node].CPD
            for example in CPD:
                CPD[example] = parameter_value

        return bn   


    def logLikelihoodIFromParams(self, variables, samples):
        # Calculates likelihood for missing data
        bn = self.BayesianNetwork

        prob_table = self.generateDataProbabilityTable()

        likelihood = decimal.Decimal('1.0')

        for sample in samples:
            sample_list = []
            missing_variable_list = []
            for k in range(len(variables)):
                variable = variables[k]
                value = sample[k]
                if (value!=""):
                    sample_list.append((variable, value))
                else:
                    missing_variable_list.append(variable)

            if len(missing_variable_list)!=0:
                #Generate variables combinations
                summ = decimal.Decimal('0.0')
                for comb in self.generateVariableValuesCombinations(missing_variable_list, bn):
                    sample_plus_missing = sample_list + comb
                    sample_frozenset = frozenset(sample_plus_missing)
                    summ = summ + prob_table[sample_frozenset]

                prob = summ
            else:
                sample_frozenset = frozenset(sample_list)
                prob = prob_table[sample_frozenset]

            likelihood = likelihood + prob.ln()
        return likelihood

    def __parallelDKESS(self, variables, data, ESS, vars_ess, bn):
        # Subrotine of Parallel version of Daphne Koller ESS computation    
        for example in data:
            # Build example with values : (variable, value)
            example_with_values = []
            for k in range(len(variables)):
                if (example[k].strip()!=""):
                    tupl = (variables[k], example[k])
                    example_with_values.append(tupl)
                    
            # Compute ESS
            for variable in vars_ess:
                for comb in ESS[variable]:
                    
                    query = list(comb)
                       
                    num = self.variableEliminationInference(query + example_with_values, [], bn)

                    if len(num)==1:
                        num = num[num.keys()[0]]
                    elif len(num)==0 :
                        num = decimal.Decimal('0.0')
                        
                    den = self.variableEliminationInference(example_with_values,[], bn)

                    if len(den)==1:
                        den = den[den.keys()[0]]
                    elif len(den)==0 :
                        den = decimal.Decimal('0.0')

                    if den!=decimal.Decimal('0.0'):
                        prob = num/den
                    else:
                        prob = decimal.Decimal('0.0')
    
                    ESS[variable][comb] = ESS[variable][comb] + prob
        return ESS
    
    

    def computeESS(self, variables, data,  comb1=None, bn=None, compute_all_prev=False, comb2=None, local=[]):
        st = time.time()
        # Compute ESS
        ESS = {}

        if compute_all_prev:
            ESS, comb1 = self.computeAllESS(variables, data, bn,  comb1, use_prob_table=False, local=local)
            ESS, comb2 = self.calculateESStoBN(ESS, data.shape[0], bn, comb2, None, nodes_to_update=local)

        else:
            # Daphne koller's EM algorithm
            #Initialize variables
            vars_ess = variables
            if len(local)>0:
                vars_ess = local
            
            ESS_structure = {}
            for variable in vars_ess:
                parents = bn[variable].Parents
                for comb in  self.generateVariableValuesCombinations([variable] + parents, bn):
                    if variable not in ESS_structure:
                        ESS_structure[variable] = {}
                    f_comb = frozenset(comb)
                    ESS_structure[variable][f_comb] = decimal.Decimal('0.0')
                           
            bn_backup = self.BayesianNetwork
            self.BayesianNetwork = None        
                        
            size = data.shape[0]        
            n = int(size/self.cores)
            jobs =[]
            
            bn = self.getDictRepresentation(bn)
                        
            job_server.set_ncpus(self.cores)
            
            for i in xrange(0, size, n):            
                chunk = data[i:i+n]
                jobs.append(job_server.submit(self.__parallelDKESS, (variables, chunk, ESS_structure, vars_ess, bn), modules=('decimal',)))
                
            job_server.wait()
            
            ESS = ESS_structure
            for job in jobs:
                #sys.stdout = open(os.devnull, "w")
                ESS_part = job() # Causes a warning
                #sys.stdout = sys.__stdout__
                for variable in ESS_part:                                       
                    for f_comb in ESS_part[variable]: 
                        ESS[variable][f_comb] += ESS_part[variable][f_comb]
                
            if self.showTime:    
                job_server.print_stats()
            #job_server.destroy()
            
            self.BayesianNetwork = bn_backup
            
        if self.showTime:
            print "Compute ESS Time: " + str(time.time() - st)
            print ""  
        return ESS, comb1, comb2

    def countMissing(self, variables, data):
        # Returns the percentual of registers with missing values in the data
        size = decimal.Decimal(data.shape[0])
        missing = decimal.Decimal('0.0')
        for example in data:
            for val in example:
                if (val.strip()==""):
                    missing+=decimal.Decimal('1.0')
                    break
        return missing/size

    def __parallelComputeAllESS(self,variables, data, combinations, prob_table, bn):
        ESS = {}
        for example in data:
            # Build example with values : (variable, value)
            example_with_values = []
            for k in range(len(variables)):
                if (example[k].strip()!=""):
                    tupl = (variables[k], example[k])
                    example_with_values.append(tupl)
        
            # For each combination of nodes
            for comb in combinations:
                comb = list(comb)
                if not(prob_table is None):
                    prob = self.inferenceFromProbTable(prob_table, comb, example_with_values, bn)
                else:                 
                    prob = self.variableEliminationInference(comb, example_with_values, bn)
        
                if len(prob)==1:
                    prob = prob[prob.keys()[0]]
                elif len(prob)==0 :
                    prob = decimal.Decimal('0.0')
                
                # Summ to the combination ESS
                f_comb = frozenset(comb)

                if (f_comb in ESS):
                    ESS[f_comb] = ESS[f_comb] + prob
                else:
                    ESS[f_comb] = prob           
        return ESS
    
    def getDictRepresentation(self, bn):
        dict_repr = {}
        for node in bn:
            dict_repr[node] = dict(bn[node])
        return dict_repr
    

    def computeAllESS(self, variables, data,  bn=None, combinations=None,use_prob_table=False, local=[]):
        # Compute All ESS
        if bn is None:
            bn = self.BayesianNetwork
                    
        bn_backup = self.BayesianNetwork 
        self.BayesianNetwork = None # workaround for pickle problem

        if use_prob_table:
            prob_table = self.generateDataProbabilityTable(bn)
        else:
            prob_table = None  
            
            
        # From variable local, include also parents of all variables in local. To compute ESS of a node, in the CPD, there are also the parents.
        local = copy.deepcopy(local)
        for node in local:
            for parent in bn[node].Parents:
                if parent not in local:
                    local.append(parent)
                           
        # Get indexes of choosen variables to compute ESS
        local_index = []    
        if len(local)!=0:
            for k in range(len(variables)):
                var = variables[k]
                if var in local:
                    local_index.append(k)
            local_data = data[:, tuple(local_index)]
            local_variables = np.array(variables)[:, tuple(local_index)].tolist()
        else:
            local_data = data
            local_variables = variables

        if combinations is None:
            #combinations = self.generateVariableValuesCombinations(bn.keys())
            combinations = set([])
            buffer = []
            count = 0
            for example in local_data:
                count+=1
                # Build example with values : (variable, value)
                example_with_values = []
                missing = []
                missing_index = []
                for k in range(len(local_variables)):
                    if (example[k].strip()==""):
                        missing.append(local_variables[k])
                        missing_index.append(k)
                    tupl = (local_variables[k], example[k])
                    example_with_values.append(tupl)
                    
                if (len(missing) > 0):                        
                    for comb in self.generateVariableValuesCombinations(missing, bn):
                        k = 0
                        for tpl in comb:
                            index = missing_index[k]
                            example_with_values[index] = tpl
                            k+=1
                        #Check before appending
                        el = tuple(example_with_values)
                        buffer.append(el)
                        if count % 100 == 0:
                            combinations = combinations.union(set(buffer))
                            buffer = []
                else:
                    # Check before appending
                    el = tuple(example_with_values)
                    buffer.append(el)
                    if count % 100 == 0:
                        combinations = combinations.union(set(buffer))
                        buffer = []

            combinations = combinations.union(set(buffer))       
        
        # For each example
        size = data.shape[0]


        n = int(size/self.cores)
        jobs =[]
        ESS = {}
        bn = self.getDictRepresentation(bn)

        job_server.set_ncpus(self.cores)
        
        for i in xrange(0, size, n):            
            chunk = data[i:i+n]
            jobs.append(job_server.submit(self.__parallelComputeAllESS, (variables, chunk, combinations, prob_table, bn), modules=('decimal',)))
        
        job_server.wait()
        for job in jobs:
            #sys.stdout = open(os.devnull, "w")
            ESS_part = job() # Causes a warning
            #sys.stdout = sys.__stdout__
            for f_comb in ESS_part: 
                if f_comb in ESS:
                    ESS[f_comb] += ESS_part[f_comb]
                else:
                    ESS[f_comb] = ESS_part[f_comb]  
        
        if self.showTime:    
            job_server.print_stats()
        #job_server.destroy()
        
        self.BayesianNetwork = bn_backup
        
        return ESS, combinations
        
    def __parallelCalculateESStoBN(self, nodes_to_update, combinations, ESS_all, bn):
        ESS={}
        for node in nodes_to_update:
            ESS[node] = {}
            parents = bn[node]['Parents']

            if combinations is None:
                combinations = {}
                combinations[node] = []
                for comb in self.generateVariableValuesCombinations(parents+[node], bn):
                    combinations[node].append(comb)
            else:
                if node not in combinations:
                    combinations[node] = []
                    for comb in self.generateVariableValuesCombinations(parents+[node], bn):
                        combinations[node].append(comb)

            for comb in combinations[node]:
                cpd = self.eliminateComplementaryLines(comb, [ESS_all], bn)[0]
                prob = sum(cpd.values())
                
                f_comb = frozenset(comb)
                ESS[node][f_comb] = prob
        


        return ESS, combinations
        

    def calculateESStoBN(self, ESS_all, len_data, bn=None, combinations=None, ESS=None,  nodes_to_update=[]):
        # Calculates the ESS for all the combinations of nodes and parents in a bn, like cpds
        if bn is None:
            bn = self.BayesianNetwork
        if ESS is None:
            ESS  = {}

        p = ProgressBar()
        count=decimal.Decimal('0.0')

        if len(nodes_to_update) ==0:
            nodes_to_update = bn.keys()
        size = len(nodes_to_update)

        n = int(size/self.cores)
        bn = self.getDictRepresentation(bn)
        
        jobs = []
        job_server.set_ncpus(self.cores)
        if n != 0:
            for i in xrange(0, size, n):
                chunk = nodes_to_update[i:i+n]
                jobs.append(job_server.submit(self.__parallelCalculateESStoBN, (chunk, combinations, ESS_all, bn), modules=('decimal',)))
        else:
            jobs.append(job_server.submit(self.__parallelCalculateESStoBN, (nodes_to_update, combinations, ESS_all, bn), modules=('decimal',)))            
        job_server.wait()
        for job in jobs:
            #sys.stdout = open(os.devnull, "w")
            ESS_part, comb_part = job() # Causes a warning
            #sys.stdout = sys.__stdout__
            
            for node in ESS_part:
                
                ESS[node] = {}
                
                if combinations is None:
                    combinations={}
                
                if node not in combinations:
                    combinations[node] = []
                    for comb in comb_part[node]:
                        combinations[node].append(comb)
                     
                for f_comb in ESS_part[node]: 
                    ESS[node][f_comb] = ESS_part[node][f_comb]
                
            del ESS_part
            del comb_part
             
        if self.showTime:
            job_server.print_stats()
        
        #job_server.destroy()
        
        del jobs
        del bn
        
        gc.collect()
        
        return ESS, combinations

    def expectedLogLikelihood(self, variables, data, ESS=None, bn=None, combinations=None):
        # Computes llikelihood when there is missing data
        if bn is None:
            bn = self.BayesianNetwork

        eLL=None
        if ESS is None:
            ESS, comb1, comb2, eLL = self.computeESS(variables, data, combinations, bn, True, None)
        
        if eLL!=None:
            return eLL
        else:    
            summ = 0
            for node in bn:
                for comb in ESS[node]:
                    parents = []
                    for tupl in comb:
                        if tupl[0] != node:
                            parents.append(tupl)
                        else:
                            node_value = tupl
                            
                    parameter = self.getProb(node_value, parents, bn[node].CPD)
                    if parameter == decimal.Decimal('0.0'): 
                        value = decimal.Decimal('0.0')
                    else:
                        value = ESS[node][comb] + parameter.ln()
                    summ += value
            return summ
        
    def expectationMaximization(self, variables, data, min_likelihood_delta, bn=None, param_ini=None, max_iter=10, local=[], computeAllESSprev=True):

        bn=self.copy(bn)
        # Expectation Maximization Algorithm
        # Random initiallization option
        if param_ini== "random":
            bn = self.randomParameterBNInitialization(bn)

        # Fixed probability param initialization
        if isinstance(param_ini, float) or isinstance(param_ini, decimal.Decimal):
            bn = self.fixedParameterBNInitialization(param_ini, bn)

        first = True
        delta = min_likelihood_delta + 1
        combinations = None
        comb2 = None
        iter = 1
        while (delta > min_likelihood_delta) and (iter <= max_iter):
            if first:
                # Expectation step
                ESS, combinations, comb2 = self.computeESS(variables, data, combinations, bn, computeAllESSprev, comb2, local=local)
                eLL = self.logLikelihood(variables, data, bn)
                first=False
            
            # Maximization step
            for node in ESS:
                parents_total = {}
                parents_ESS = self.summout(ESS[node], node, bn)
                for comb in ESS[node]:
                    parents = []
                    for tupl in comb:
                        if tupl[0]!=node:
                            parents.append(tupl)

                    f_parents = frozenset(parents)
                    # Laplace Smoothing                   
                    parameter = ((len(bn[node].Values) * ESS[node][comb]) + self.smoothing_param)/((len(bn[node].Values) * parents_ESS[f_parents]) + (self.smoothing_param * decimal.Decimal(len(bn[node].Values))))             
                                                
                    bn[node].CPD[comb] = parameter
                        
            last = eLL
            # Expectation step
            ESS, combinations, comb2 = self.computeESS(variables, data, combinations, bn, computeAllESSprev, comb2, local=local)
            eLL = self.logLikelihood(variables, data, bn)
            
            
            # Delta calculation
            delta = eLL - last

            print "Likelihood increase: " + str(delta)
            iter += 1

        return ESS, bn

    def __is_bottom_marked(self, nodes_status, node):
        # Checks if the bottom of a node is marked in bayes ball
        bottom_marked = False
        if ((nodes_status[node]==2) or (nodes_status[node]==3)):
            bottom_marked = True
        return bottom_marked

    def __is_top_marked(self, nodes_status, node):
        # Checks if the top of a node is marked in bayes ball
        top_marked = False
        if ((nodes_status[node]==1) or (nodes_status[node]==3)):
            top_marked = True
        return top_marked

    def bayesBall(self, query_nodes, evidence_nodes, bn=None):
        # Implementation of the BayesBall algorithm
        # The query node should not be in the evidence
        if bn is None:
            bn = self.BayesianNetwork
        # Possible Status:
        # 0 - Not Visited
        # 1 - Marked on the top
        # 2 - Marked on the bottom
        # 3 - Marked on the top and on the bottom

        # Set not visited status to every node
        nodes_status = {}
        for node in bn:
            nodes_status[node] = 0

        nodes_to_be_visited = []
        nodes_already_visited = []

        # From where to visited
        # True - From Parent
        # False - From Child
        from_where_to_visit =[]
        # Sets each node to be visited from child
        for node in query_nodes:
            nodes_to_be_visited.append(node)
            from_where_to_visit.append(False)

        # While there are nodes to be visited
        while len(nodes_to_be_visited)>0:
            choosen_node = nodes_to_be_visited.pop()
            from_where_to_visit_choosen_node = from_where_to_visit.pop()
            choosen_node_parents = bn[choosen_node].Parents
            #Get children node from the choosen node
            choosen_node_children = []
            for node in bn:
                if choosen_node in bn[node].Parents:
                    choosen_node_children.append(node)
            # Mark as visited
            nodes_already_visited.append(choosen_node)
            
            if (choosen_node not in evidence_nodes) and (from_where_to_visit_choosen_node == False):
                # if the top is not marked
                if not (self.__is_top_marked(nodes_status,choosen_node)):
                    #Mark top
                    if nodes_status[choosen_node] == 2:
                        nodes_status[choosen_node] = 3
                    else:
                        nodes_status[choosen_node] = 1

                    #Schedule parents to be visited
                    for node in choosen_node_parents:
                        nodes_to_be_visited.append(node)
                        from_where_to_visit.append(False)

                 # if the bottom is not marked
                if  not (self.__is_bottom_marked(nodes_status,choosen_node)):
                    #Mark bottom
                    if nodes_status[choosen_node] == 1:
                        nodes_status[choosen_node] = 3
                    else:
                        nodes_status[choosen_node] = 2

                    #Schedule children to be visited
                    for node in choosen_node_children:
                        nodes_to_be_visited.append(node)
                        from_where_to_visit.append(True)

            # if the visit is from a parent
            if (from_where_to_visit_choosen_node == True):
                # If choosen_node is  in the evidence nodes and its top is not marked
                if ((choosen_node in evidence_nodes) and (not (self.__is_top_marked(nodes_status,choosen_node)))):
                    #Mark top
                    if nodes_status[choosen_node] == 2:
                        nodes_status[choosen_node] = 3
                    else:
                        nodes_status[choosen_node] = 1

                    #Schedule parents to be visited
                    for node in choosen_node_parents:
                        nodes_to_be_visited.append(node)
                        from_where_to_visit.append(False)

                # if choosen_node not in evidence and its bottom is not marked
                if (not(choosen_node in evidence_nodes) and (not (self.__is_bottom_marked(nodes_status,choosen_node)))):
                    #Mark bottom
                    if nodes_status[choosen_node] == 1:
                        nodes_status[choosen_node] = 3
                    else:
                        nodes_status[choosen_node] = 2

                    #Schedule children to be visited
                    for node in choosen_node_children:
                        nodes_to_be_visited.append(node)
                        from_where_to_visit.append(True)

        # Get results
        irrelevant_nodes = []
        requisite_probability_nodes = []
        requisite_observation_nodes = list(set(nodes_already_visited).intersection(set(evidence_nodes)))

        for node in nodes_status:
            # Irrelevant nodes - not marked on the bottom
            if not (self.__is_bottom_marked(nodes_status, node)):
                irrelevant_nodes.append(node)

            # Requisite probability nodes - marked on the top
            if self.__is_top_marked(nodes_status, node):
                requisite_probability_nodes.append(node)

        
        return irrelevant_nodes, requisite_probability_nodes, requisite_observation_nodes, list(set(nodes_already_visited))

    def existEdge(self, nodeA, nodeB, bn):
        # Check if an edge between nodeA and nodeB exists
        if bn is None:
            bn = self.BayesianNetwork
        exist = False
        if ((nodeB in bn[nodeA].Parents) or (nodeA in bn[nodeB].Parents)):
            exist = True
        return exist

    def removeEdge(self, nodeA, nodeB, bn):
        # Remove the edge of the bayesian network structure that connects nodeA and nodeB
        if bn is None:
            bn = self.BayesianNetwork

        if (nodeA in bn[nodeB].Parents):
            # Remove edge
            bn[nodeB].Parents.remove(nodeA)
        elif (nodeB in bn[nodeA].Parents):
            # Remove edge
            bn[nodeA].Parents.remove(nodeB)

        return bn


    def __removeNode(self, choosen_node, bn_structure):
        # Removes a node from a graph
        del bn_structure[choosen_node]
        for node in bn_structure:
            if choosen_node in bn_structure[node]:
                bn_structure[node].remove(choosen_node)
        return bn_structure

    def __getLeaves(self, bn_structure):
        # Get the leaves of a graph
        # Nodes that don't have children are nodes that are not parents of any other node
        nodes = bn_structure.keys()
        parents = []
        for node in nodes:
            parents  = parents + bn_structure[node]
        return list(set(nodes).difference(set(parents)))

    def checkAciclicity(self, bn):
        # Check if the graph is still aciclic.
        if bn is None:
            bn = self.BayesianNetwork
        # Get Bayesian Network structure
        bn_structure = {}
        for node in bn:
            bn_structure[node] = bn[node].Parents

        leaves = []
        aciclic = True
        while (True):
            # Check if exist nodes
            if len(bn_structure)==0:
                break

            # Get the leaves of the graph
            leaves =  self.__getLeaves(bn_structure)

            # Check if exist leaves
            if len(leaves)==0:
                aciclic = False
                break

            # Remove leaf
            leaf = leaves.pop()
            bn_structure = self.__removeNode(leaf, bn_structure)

        return aciclic

    def addEdge(self, nodeA, nodeB, bn):
        #Add an edge from nodeA to nodeB
        if bn is None:
            bn = self.BayesianNetwork
        bn[nodeB].Parents.append(nodeA)
        return bn


    def invertEdge(self, nodeA, nodeB, bn):
        # Invert the direction of an edge from nodeA to nodeB
        if bn is None:
            bn = self.BayesianNetwork

        if (nodeA in bn[nodeB].Parents):
            bn[nodeB].Parents.remove(nodeA)
            bn[nodeA].Parents.append(nodeB)
        elif (nodeB in bn[nodeA].Parents):
            bn[nodeA].Parents.remove(nodeB)
            bn[nodeB].Parents.append(nodeA)

        return bn

    def generateVariablesCombinations(self, variables_list, list_size):
        comb_list = copy.deepcopy(variables_list)
        for i in range(2, list_size+1):
            # Combine 2-2
            end = len(comb_list)
            for k in range(end):
                item_comb = comb_list.pop(0)
                if isinstance(item_comb, str) or isinstance(item_comb, unicode):
                    item_comb = [item_comb]
                last = item_comb[-1]
                idx = variables_list.index(last)
                for j in range(idx+1, len(variables_list)):
                    item_var = [variables_list[j]]
                    tmp = copy.deepcopy(item_comb)
                    tmp = tmp + item_var
                    comb_list.append(tmp)
        return comb_list

    def copy(self, bn, new_bn=None, cpds_only=False, nodes_to_copy=[]):
        # Copy a bayesian network. Copy part of a bayesian network to another one.
        if new_bn is None:
            new_bn = Storage()
        
        if len(nodes_to_copy)==0:
            nodes_to_copy = bn.keys()
        
        for node in nodes_to_copy:
            if not cpds_only:
                new_bn[node] = Storage()
                new_bn[node].Values = copy.deepcopy(bn[node].Values)
                new_bn[node].Parents = copy.deepcopy(bn[node].Parents)
            new_bn[node].CPD = copy.deepcopy(bn[node].CPD)
        return new_bn


    def generateNetNeighbours(self, nodes_to_consider=None, external_edges_to_all=None, mutual_information_filter=None, bn=None, dataset_size=None, ESS_all=None, edges_to_consider=None, max_neighbours=50):
        # Generate all neighbour structure of a given bayesian network
        # A neighbour network is one that considers only one change (add edge, remove edge or  invert egde)
        # The nodes_to_consider are the nodes which we want the edges between then to be changed.
        # external_edges_to_all is a list of nodes from where it will be created neighbours that includes conections to all nodes
        # edges_to_consider allows generating neighbours by only changing a group of edges (associations between two nodes)

        # Get the structure of the network we have
        bn_structure = {}
        for node in bn:
            bn_structure[node] = Storage()
            bn_structure[node].Parents = bn[node].Parents
            bn_structure[node].Values = bn[node].Values

        # Generate possible node combinations
        if nodes_to_consider is None:
            nodes = bn_structure.keys()
        else:
            nodes = nodes_to_consider

        # If user has determined the edges to consider, use his edges
        if edges_to_consider == None:
            tmp_combs = self.generateVariablesCombinations(nodes,2)
        else:
            tmp_combs = edges_to_consider

        # Add more possible combinations (considering other nodes in the network) (External edges)
        if external_edges_to_all != None:
            print "       Adding external edges..."
            external = set(bn.keys()).difference(set(nodes))
            for node1 in external_edges_to_all:
                for node2 in external:
                    tmp_combs.append([node1, node2]);
                    
        # Mutual information filter
        if mutual_information_filter!=None:
            print "       Applying mutual information filter..."
            combs = []
            mi_combs = {}
            for comb in tmp_combs:
                mi = self.expectedMutualInformation(comb[0], comb[1], dataset_size, ESS_all, bn)
                combs.append(tuple(comb))
                mi_combs[tuple(comb)] = mi
            ordered_by_mi = sorted(combs, key=lambda cb: mi_combs[cb])
            size = int((len(ordered_by_mi) * mutual_information_filter).to_integral_exact(rounding=decimal.ROUND_CEILING))
            combs = ordered_by_mi[:size]
        else:
            combs = tmp_combs

        # Shuffle list (everyday I'm shuff-ling-ling)
        random.shuffle(combs)
        
        # Generate neighbours
        print "       Generating neighbours..."
        counter = 0
        for comb in combs:          
            # Evaluate involved nodes situation
            node1 = comb[0]
            node2 = comb[1]
            new_bn_structure1 = None
            new_bn_structure2 = None
            # Check if exist an edge
            if self.existEdge(node1, node2, bn_structure):
                # first possible neighbour
                new_bn_structure1 = self.copy(bn_structure)
                new_bn_structure1 = self.invertEdge(node1, node2, new_bn_structure1)
                # check if is a valid structure
                if (self.checkAciclicity(new_bn_structure1)):
                    neighbour_structure = new_bn_structure1
                    affected_node = [node1, node2]
                    print "\tNeighbour change: Inverted (" + node1  + "," + node2 + ")"
                    yield neighbour_structure, affected_node
                    counter+=1
                    if counter >= max_neighbours:
                        break
                # second possible neighbour
                new_bn_structure2 =  self.copy(bn_structure)
                new_bn_structure2 = self.removeEdge(node1, node2, new_bn_structure2)
                # check if is a valid structure
                if (self.checkAciclicity(new_bn_structure2)):
                    neighbour_structure = new_bn_structure2
                    affected_node = [node1, node2]
                    print "\tNeighbour change: Removed (" + node1  + "," + node2 + ")"
                    yield neighbour_structure, affected_node
                    counter+=1
                    if counter >= max_neighbours:
                        break
            else:
                # first possible neighbour
                new_bn_structure1 =  self.copy(bn_structure)
                new_bn_structure1 = self.addEdge(node1, node2, new_bn_structure1)
                # check if is a valid structure
                if (self.checkAciclicity(new_bn_structure1)):
                    neighbour_structure = new_bn_structure1
                    affected_node = [node2]
                    print "\tNeighbour change: Added (" + node1  + "," + node2 + ")"
                    yield neighbour_structure, affected_node
                    counter+=1
                    if counter >= max_neighbours:
                        break
                
                # second possible neighbour
                new_bn_structure2 = self.copy(bn_structure)
                new_bn_structure2 = self.addEdge(node2, node1, new_bn_structure2)
                if (self.checkAciclicity(new_bn_structure2)):
                    neighbour_structure = new_bn_structure2
                    affected_node = [node1]
                    print "\tNeighbour change: Added (" + node2  + "," + node1 + ")"
                    yield neighbour_structure, affected_node
                    counter+=1
                    if counter >= max_neighbours:
                        break
                    
                    
    def createVariable(self):
        # Creates a variable to be appended to a network
        var = Storage()
        var.CPD = {}
        var.Values = []
        var.Parents = []
        return var
    
    def generateDAHVIneighbours(self, revision_points, values, bn):
        # Generate the candidate networks to be evaluated by DAHVI algorithm
        for rev_point in revision_points:
            neighbour_network = self.copy(bn)
            rev_parents = bn[rev_point].Parents
            
            # Give a name to the hidden variable
            name = "Hidden"
            counter = 1
            while ((name+str(counter)) in bn.keys()):
                counter+=1
            
            # Generate Hidden variable
            hidden_var = self.createVariable()
            hidden_var_name = name+str(counter)
            neighbour_network[hidden_var_name] = hidden_var
            neighbour_network[hidden_var_name].Values = values
            neighbour_network[hidden_var_name].Parents = rev_parents
            neighbour_network[rev_point].Parents = [hidden_var_name]
            
            # Initialize CPD of the hidden variables with zero
            neighbour_network = self.randomParameterBNInitialization(neighbour_network, [hidden_var_name, rev_point], overwrite=True)
            
            yield neighbour_network
                
            

    def inferenceFromData(self, variables, data, query, evidence, bn=None, prob_table=None):
        # Do probabilistic inference by counting from data
        if bn is None:
            bn = self.BayesianNetwork

        # Generate probability table: all possible combinations of values, and its probabilities
        if prob_table is None:
            prob_table = self.generateDataProbTableFromDataFreq(variables, data, bn)

        # Format the query
        query_list, ev_list = self.buildQuery(query, evidence)

        # Return the variables in query and evidence
        vars_with_values  = query_list + ev_list
        variables_in_query = []
        for var_tpl in vars_with_values:
            variables_in_query.append(var_tpl[0])

        # Filter results
        tmp = []
        tmp.append(prob_table)
        prob_table = self.eliminateComplementaryLines(evidence, tmp, bn);
        total = decimal.Decimal('0.0')
        for prob_key in prob_table[0]:
            total  += prob_table[0][prob_key]

        # Sum over variables that are not included
        to_summout =  list((set(variables)).difference(set(variables_in_query)))
        cpd = prob_table[0]
        for variable in to_summout:
            cpd = self.summout(cpd, variable, bn)

        not_elim_list = []
        for query_tpl in query_list:
            if query_tpl[1] != self.any_value:
                not_elim_list.append(query_tpl)

        if len(not_elim_list) != 0:
            tmp = []
            tmp.append(cpd)
            cpd = self.eliminateComplementaryLines(not_elim_list, tmp, bn)

        num = decimal.Decimal('0.0')
        for prob_key in cpd[0]:
            num  += cpd[0][prob_key]

        return num/total

    def __parallelLocalLogLikelihood(self, variables, data, variable, parents, bn=None):
        # Parallel part of computing decomposable log likelihood
        local_likelihood = decimal.Decimal("0.0")

        for example in data:
            query = []
            evidence = []
            for k in range(len(example)):
                if example[k].strip() != "":
                    if variables[k]==variable:                   
                        query.append((variables[k], example[k]))
                    if variables[k] in parents:
                        evidence.append((variables[k], example[k]))
                           
            inference = self.variableEliminationInference(query, evidence, bn)
            if len(inference)==1:
                inference = inference[inference.keys()[0]]
            elif len(inference)==0 :
                inference = decimal.Decimal('0.0')
              
            local_likelihood += inference.ln()
        
        return local_likelihood


    def localLogLikelihood(self, variables, data, variable, bn):
        # Computes the likelihood score for a variable and its parents, assuming it is decomposable (complete data case)
        size = data.shape[0]
        n = int(size/self.cores)

        parents = bn[variable].Parents

        bn = self.getDictRepresentation(bn)
        
        bn_backup = self.BayesianNetwork 
        self.BayesianNetwork = None # workaround for pickle problem

        jobs = []
        job_server.set_ncpus(self.cores)
        if n != 0:
            for i in xrange(0, size, n):
                chunk = data[i:i+n]
                jobs.append(job_server.submit(self.__parallelLocalLogLikelihood, (variables, chunk, variable, parents, bn), modules=('decimal',)))
        else:
            jobs.append(job_server.submit(self.__parallelLocalLogLikelihood, (variables, data, variable, parents, bn), modules=('decimal',)))            
        job_server.wait()
        
        total_likelihood = 0
        for job in jobs:
            #sys.stdout = open(os.devnull, "w")
            likelihood_part = job() # Causes a warning
            #sys.stdout = sys.__stdout__
            total_likelihood+=likelihood_part
            
                         
        if self.showTime:
            job_server.print_stats()
        
        self.BayesianNetwork = bn_backup

        return total_likelihood


    def logLikelihoodDecomposable(self, variables, data, score_cache, bn, update_variables=[]):
        # Computes loglikelihood score in the case of a complete dataset
        if (update_variables==[]):
            update_variables = bn.keys()
            
        if score_cache is None:
            score_cache = {}
            
        loglikelihood = decimal.Decimal("0.0")
        for variable in update_variables:
            score_cache[variable] = self.localLogLikelihood(variables, data, variable, bn)
            loglikelihood += score_cache[variable]
            
        return loglikelihood, score_cache
    
    
    def __parallelLogLikelihood(self, variables, data, bn=None):
        
        loglikelihood = 0
        for example in data:
            query = []
            for k in range(len(example)):
                if example[k].strip() != "":
                    query.append((variables[k], example[k]))
                    
            inference = self.variableEliminationInference(query, [], bn)
            if len(inference)==1:
                inference = inference[inference.keys()[0]]
            elif len(inference)==0 :
                inference = decimal.Decimal('0.0')
   
            loglikelihood += inference.ln()

        return loglikelihood

    
    def logLikelihood(self, variables, data, bn=None):
        # Computes loglikelihood score        
        size = data.shape[0]
        n = int(size/self.cores)

        bn = self.getDictRepresentation(bn)
        
        bn_backup = self.BayesianNetwork 
        self.BayesianNetwork = None # workaround for pickle problem

        jobs = []
        job_server.set_ncpus(self.cores)
        if n != 0:
            for i in xrange(0, size, n):
                chunk = data[i:i+n]
                jobs.append(job_server.submit(self.__parallelLogLikelihood, (variables, chunk, bn), modules=('decimal',)))
        else:
            jobs.append(job_server.submit(self.__parallelLogLikelihood, (variables, data, bn), modules=('decimal',)))            
        job_server.wait()
        
        total_likelihood = 0
        for job in jobs:
            #sys.stdout = open(os.devnull, "w")
            likelihood_part = job() # Causes a warning
            #sys.stdout = sys.__stdout__
            total_likelihood+=likelihood_part
            
                         
        if self.showTime:
            job_server.print_stats()
        
        self.BayesianNetwork = bn_backup
        return total_likelihood
    
    

    def generateRandomStructure(self, variables, prob=decimal.Decimal('0.5'), max_parents=None):
        # Generate a random bayesian structure from the given variables, whith probability prob of existing an edge
        bn = Storage()
        for variable in variables:
            if variable not in bn:
                bn[variable] = Storage()
                bn[variable].Parents = []
                bn[variable].Values = self.BayesianNetwork[variable].Values

        for variable1 in variables:
            for variable2 in variables:
                exist_edge = self.simpleSampling([True, False], [prob, decimal.Decimal('1.0')-prob])
                if (exist_edge):
                    is_parent = self.simpleSampling([True, False], [prob, decimal.Decimal('1.0')-prob])
                    if (is_parent):
                        if (not (variable1 in bn[variable2].Parents)) and (len(bn[variable2].Parents) < max_parents):
                            bn[variable2].Parents.append(variable1)
                            if not (self.checkAciclicity(bn)):
                                bn[variable2].Parents.remove(variable1)
                    else:
                        if (not (variable2 in bn[variable1].Parents)) and (len(bn[variable1].Parents) < max_parents):
                            bn[variable1].Parents.append(variable2)
                            if not (self.checkAciclicity(bn)):
                                bn[variable1].Parents.remove(variable2)
        return bn
    
    def generateStructureWithoutEdges(self, variables):
        # Generate a structure without edges, only the nodes, as if all the variables were independent
        bn = self.generateRandomStructure(variables, decimal.Decimal("0.0"), None)
        return bn
        

    def greedyHillClimbing(self, variables, data, complete_data=True, score_type="likelihood", ESS_all=None, initial_structure=None, edges_to_consider=None, max_iter=decimal.Decimal("inf"), class_variable=None, nodes_to_consider=None, external_edges_to_all=None, max_neighbours=50, path_partial_out=None):
        # Greedy Hill Climbing, in its simpler form, for searching structures.
        # Data can be complete or incomplete
        # Score can be likelihood for loglikelihood score, bic for BIC score, cll for conditional loglikelihood
        if initial_structure is None:
            bn = self.BayesianNetwork
        else:
            bn = self.copy(initial_structure)
        
        # Calculate initial score
        if score_type =="likelihood":
            score_now = self.logLikelihood(variables, data, bn)
        if score_type =="bic":
            score_now = self.BICScore(variables, data, complete_data, bn, None)
        if score_type =="cll":
            score_now = self.conditionalLogLikelihood(variables, data, class_variable, bn)
        if score_type =="correct_class":
            score_now = self.countCorrectClassifications(variables, data, class_variable, decimal.Decimal("0.5"), bn)[0]
        if score_type == "aic":
            score_now = self.AIC(variables, data, bn, "cll", class_variable)
                
        score_before = decimal.Decimal("-inf")
        # while score is still increasing
        progress = True
        best_bn = bn
        #best_bn_ESS = ESS_ini
        iter = 0
        while (progress and (iter < max_iter)):
            iter += 1
            progress = False
            score_before = score_now
            bn = best_bn
            #ESS_ini = copy.deepcopy(best_bn_ESS)
                                 
            # For each neighbour of the bn
            neighbours_count = 0
            criteria = False
            for ret  in self.generateNetNeighbours(nodes_to_consider, external_edges_to_all, None, bn, None, None, edges_to_consider, max_neighbours=max_neighbours):  
                st = time.time()  
                neighbours_count += 1
                neighbour_structure = ret[0]
                affected_node = ret[1]
                
                if complete_data:
                    if score_type =="likelihood":
                        score = self.logLikelihood(variables, data, neighbour_structure)
                else:
                    neighbour_structure = self.copy(bn, neighbour_structure, True)
                    neighbour_structure = self.fixedParameterBNInitialization(decimal.Decimal('0.0'), neighbour_structure, affected_node, True)
                    neighbour_structure = self.estimateBNParamsFromESS(ESS_all, neighbour_structure,  affected_node)
                    
                    if score_type =="likelihood":
                        score = self.logLikelihood(variables, data, neighbour_structure)
                    if score_type =="bic":
                        score = self.BICScore(variables, data, complete_data, neighbour_structure, None)
                    if score_type=="cll":
                        score = self.conditionalLogLikelihood(variables, data, class_variable, neighbour_structure)
                    if score_type=="correct_class":
                        score = self.countCorrectClassifications(variables, data, class_variable, decimal.Decimal("0.7"), neighbour_structure)[0]
                    if score_type=="aic":
                        score = self.AIC(variables, data, neighbour_structure, "cll", class_variable)

                print "\t{0}) Score: {1} Time: {2}".format(str(neighbours_count), str(score), str(time.time()-st))
                
                if score_type == "aic":
                    criteria = (score < score_now)
                else:
                    criteria = (score > score_now)
                
                if (criteria):
                    print "\t Found a better score!"
                    score_now = score
                    best_bn = neighbour_structure
                    progress = True
                    
                #del ESS
                del neighbour_structure
                del affected_node
            
            if path_partial_out != None:
                self.export_encog(path_partial_out, best_bn)
                
            print "Neighbours evaluated: " + str(neighbours_count)
            print "Structure score increase: " +  str(score_now - score_before) + "; from: " + str(score_before) + " to: " + str(score_now)

        return best_bn

    def generateBNStructureVisualization(self, path_output, bn=None):
        # Generate a picture with a visualization of the network
        if bn is None:
            bn = self.BayesianNetwork

        graph = pydot.Dot(graph_type='digraph')

        variables = bn.keys()
        for variable in variables:
            node =  pydot.Node(variable, style="filled", fillcolor="white")
            graph.add_node(node)
            for parent in bn[variable].Parents:
                graph.add_edge(pydot.Edge(parent, node))

        graph.write_png(path_output)
        
        
    def numberOfParameters(self, bn):
        # Counts the number of parameters in the bayesian network
        sum = 0
        prod = 0
        for node in bn:
            prod = len(bn[node].Values)
            for parent in bn[node].Parents:
                prod = prod * len(bn[parent].Values)
            sum += prod
            
        return sum
    

    def BICScore(self, variables, data, complete_data=True,  bn_structure=None, ESS=None):
        # calculates the BIC score, for complete or incomplete data
        # Compute dimension of the graph, number of independent parameters
        if bn_structure is None:
            bn_structure = self.BayesianNetwork

        dim_g =  self.numberOfParameters(bn_structure)
        bic = self.logLikelihood(variables, data, bn_structure)

        m = decimal.Decimal(data.shape[0])
        bic = bic - (m.ln()/2) * dim_g

        return bic

    def inferenceFromProbTable(self, prob_table, query, evidence, bn=None):
        if bn is None:
            bn  = self.BayesianNetwork

        # Calculates the probabilities given that every possibility, every mix of values, is already calculated in a probability table
        query_var_values =  {}
        ev_var_values = {}

        query, evidence = self.buildQuery(query, evidence)

        for tpl  in query:
            query_var_values[tpl[0]] = tpl[1]

        for tpl in evidence:
            ev_var_values[tpl[0]] = tpl[1]

        # Check if query and evidence are oposed
        for var in query_var_values:
            if var in ev_var_values:
                if (query_var_values[var] != ev_var_values[var]):
                    if (query_var_values != self.any_value):
                        return {}

        vars = bn.keys()
        vars = set(vars).difference(set(query_var_values.keys() + ev_var_values.keys()))

        if len(evidence)==0:
            cpd =  self.eliminateComplementaryLines(query+evidence, [prob_table], bn)[0]
        else:   
            cpd =  self.eliminateComplementaryLines(evidence, [prob_table], bn)[0]

        for var in vars:
            cpd = self.summout(cpd, var, bn)

        if len(evidence)!=0:
            numerator  = cpd
            vars = set(query_var_values.keys()).difference(ev_var_values.keys())
            for var in vars:
                cpd = self.summout(cpd, var, bn)
            denominator  = cpd
            cpd = self.normalize(numerator, denominator)

        keep_vals = []
        for tpl in (query+evidence):
            if tpl[1] != self.any_value:
                keep_vals.append(tpl)

        if len(keep_vals)!=0:
            cpd_set = self.eliminateComplementaryLines(keep_vals, [cpd], bn)
        else:
            cpd_set = [cpd]

        return cpd_set[0]

    def estimateBNParamsFromESS(self, ESS, bn_structure, nodes_to_consider=[]):
        # Estimates  BN parameters in the CPDs using the ESS. Updates the parameters in the BN CPDs
        if bn_structure is None:
            bn_structure = self.BayesianNetwork
        else:
            bn_structure = self.copy(bn_structure)
            
        for node in nodes_to_consider:
            cpd  = bn_structure[node].CPD
            for comb in cpd:
                parents = []
                for tpl in comb:
                    if tpl[0]!=node:
                        parents.append(tpl)
                            
                den_CPD = self.eliminateComplementaryLines(parents, [ESS], bn_structure)[0]       
                den = sum(den_CPD.values())
                num_CPD = self.eliminateComplementaryLines(comb, [den_CPD], bn_structure)[0]
                num = sum(num_CPD.values())
                
                #Laplace Smoothing
                bn_structure[node].CPD[comb] = ((num*len(bn_structure[node].Values)) + self.smoothing_param)/((den*len(bn_structure[node].Values)) + self.smoothing_param*len(bn_structure[node].CPD))
                
                
        return bn_structure

    def printBN(self, bn=None):
        # Print a full bayesian network, not only the CPDs, but the variables, and the values also.
        if bn is None:
            bn = self.BayesianNetwork

        for node in bn:
            print "===================================================== Node: " + node  + "(" + str(len(bn[node].Values))+ " Values) ====================================================="
            parents_arr = []
            for parent in bn[node].Parents:
                parents_arr.append(parent + "("  + str(len(bn[parent].Values)) +  " Values)")
            parents_str = ", ".join(parents_arr)
            print "Parents:" + parents_str
            values_str = ", ".join(bn[node].Values)
            print "Possible Values: " + values_str
            print "CPD: (" +  str(len(bn[node].CPD)) +" entries)"
            self.printCPD(bn[node].CPD)
            print ""


    def structuralExpectationMaximization(self, variables, data, output, min_likelihood_delta=decimal.Decimal('0.01'), structure_search="greedyHC", search_options={}, bn=None, mutual_information_filter=None, max_iter=10, score_type="bic", class_variable=None, sem_max_iter=decimal.Decimal("Infinity")):
        # Structural EM algorithm (SEM)
        if bn is None:
            bn = self.BayesianNetwork

        print "Started SEM."
        # While don't converge
        bn_structure =  self.copy(bn)
        # Initial BIC score calculation
        print "Initial score calculation..."
        if score_type == "bic":
            last = self.BICScore(variables, data, False, bn_structure)
        elif score_type == "cll":
            last = self.conditionalLogLikelihood(variables, data, class_variable, bn_structure)
        elif score_type == "aic":
            last = self.AIC(variables, data, bn_structure, "cll", class_variable)
        
        print "Score type is " + score_type +". Initial score is " + str(last)    
        delta = min_likelihood_delta + 1
        iteration = 1
        while (delta > min_likelihood_delta) and (iteration <= sem_max_iter):
            print ""
            print "Iteration: " + str(iteration)
            last_bn = self.copy(bn_structure)
            # Expectation Maximization (get params, ESS)
            if structure_search == "BBSL":
                if "max_nodes" not in search_options:
                    search_options["max_nodes"] = None
                local = list(self.BBSL_StepI(variables, data, search_options["class_variable"], search_options["classification_threshold"], 0.1,bn_structure, search_options["max_nodes"]))
                print "Bayes Ball choosen nodes: " + ",".join(local)                                         
            else:
                local=[]    
            
            print "Expectation maximization..."    
            ESS_all, bn_structure = self.expectationMaximization(variables, data, min_likelihood_delta, bn_structure, None, max_iter, local=local, computeAllESSprev=False)
            
            #ESS_all, combinations = self.computeAllESS(variables, data, bn_structure, None, False, local=local)
            #ESS_all, comb1, comb2 = self.computeESS(variables, data, None, bn_structure, False, None, local=local)
            
            # Search structure (use ESS (data) to score structures)
            print "Running structure search..."
            if structure_search=="greedyHC":
                if "path_partial_out" not in search_options:
                    search_options["path_partial_out"] = None
                bn_structure = self.greedyHillClimbing(variables, data, False, "aic", ESS_all, bn_structure, edges_to_consider=None, max_iter=10,class_variable=search_options["class_variable"], max_neighbours=decimal.Decimal("Infinity"), path_partial_out=search_options["path_partial_out"])
            elif structure_search == "BBSL":
                bn_structure= self.BBSL(variables, data, search_options["class_variable"], search_options["classification_threshold"], ESS_all, 'aic',mutual_information_filter ,  bn_structure, search_options["external_edges"],0.1, set(local))

            # Maximization (calculate parameters for the new structure from ESS)
            print "Maximization step..."
            if structure_search == "BBSL":
                nodes_to_consider = local
            else:
                nodes_to_consider=bn_structure.keys()
            bn_structure = self.estimateBNParamsFromESS(ESS_all, bn_structure, nodes_to_consider)
            

            # BIC score (evaluate)
            if score_type == "bic":
                score = self.BICScore(variables, data, False, bn_structure)
            elif score_type == "cll":
                score = self.conditionalLogLikelihood(variables, data, class_variable, bn_structure)
            elif score_type == "aic":
                score = self.AIC(variables, data, bn_structure, "cll", class_variable)
                
            # Delta calculation
            delta = score - last
            print "Score increase: " + str(delta) + ". From  " + str(last) + " to " + str(score)
            last = score
            iteration += 1

            if score_type == "aic":
                delta = -delta

            # If the new structure is worse; should only happen in the first iteration
            if (delta < 0):
                bn_structure = last_bn
                                
        final_structure = self.copy(bn_structure)

        if structure_search == "BBSL":
            print "Smoothing parameters... "
            final_structure = self.laplaceSmoothBNParams(final_structure, decimal.Decimal("0.01"))
            print "Finished smoothing."

        print "Finished SEM."
        #self.generateBNStructureVisualization(output, final_structure)
        
        return final_structure


    def compareBNStructures(self, bn_1, bn_2, variables=None, data=None, class_variable=None, classification_threshold=decimal.Decimal("0.5")):
        # Compare two bayesian network structures, both bayesian network structures need to have the same nodes
        # Present nodes
        # Number of added edges
        # Number of removed edges
        # Number of inverted edges
        
        inverted = 0
        added = 0
        removed = 0
        inverted_list = []
        added_list = []
        removed_list = []
        # Compare the networks
        nodes = bn_1.keys()
        for i in range(len(nodes)):
            node_1 = nodes[i]
            for j in range(i+1,len(nodes)):
                node_2 = nodes[j]
                if (node_1 != node_2):
                    exist_edge_bn_1 = self.existEdge(node_1, node_2, bn_1)
                    exist_edge_bn_2 = self.existEdge(node_1, node_2, bn_2)

                    if exist_edge_bn_1 and exist_edge_bn_2:
                        # Check if inverted
                        if ((node_1 in bn_1[node_2].Parents) and (node_2 in bn_2[node_1].Parents)) or ((node_2 in bn_1[node_1].Parents) and (node_1 in bn_2[node_2].Parents)):
                            inverted_list.append("("+node_1+", " + node_2+")")
                            inverted+=1
                    else:
                        if exist_edge_bn_1 and (not exist_edge_bn_2):
                            removed_list.append("("+node_1+", " + node_2+")")
                            removed+=1

                        if  (not exist_edge_bn_1) and exist_edge_bn_2:
                            added_list.append("("+node_1+", " + node_2+")")
                            added+=1
                            
        print ""
        print "=====================COMPARISON RESULT========================="

        print "Added edges: " + str(added)
        print "Added edges: " + ",".join(added_list)
        print "Removed edges: "  + str(removed)
        print "Removed edges: " + ",".join(removed_list)
        print "Inverted edges: " + str(inverted)
        print "Inverted edges: " + ",".join(inverted_list)
        print "Structural Hamming Distance: " + str(added + inverted + removed)

        if not(variables is None):        
            # Scores for incomplete data
            # Expected Likelihood score        
            likelihood_bn1 = self.logLikelihood(variables, data, bn_1)
            print "Likelihood score from first BN: " + str(likelihood_bn1)
            likelihood_bn2 = self.logLikelihood(variables, data, bn_2)
            print "Likelihood score from second BN: " + str(likelihood_bn2)
            print "Likelihood score increase: " + str(likelihood_bn2 - likelihood_bn1)
    
            bic_bn1 = self.BICScore(variables, data, False,  bn_1, ESS=None)
            print "BIC score from first BN: " + str(bic_bn1)
            bic_bn2 = self.BICScore(variables, data, False, bn_2, ESS=None)
            print "BIC score from second BN: " + str(bic_bn2)
            print "BIC score increase: " + str(bic_bn2 - bic_bn1)
                   
            if class_variable != None:
                # CLL
                cond_1  = self.conditionalLogLikelihood(variables, data, class_variable, bn_1)
                print "Conditional LogLikelihood for node " + class_variable + " in the first BN: " + str(cond_1)
                cond_2 = self.conditionalLogLikelihood(variables, data, class_variable, bn_2)
                print "Conditional LogLikelihood for node " + class_variable + " in the second BN " + str(cond_2)
                print "Conditional LogLikelihood increase for node " +  class_variable + " :"  + str(cond_2 - cond_1)
                
                #AIC
                aic_1 = self.AIC(variables, data, bn_1, "cll", class_variable)
                print "AIC score with CLL for node " + class_variable + " in the first BN: " + str(aic_1)
                aic_2 = self.AIC(variables, data, bn_2, "cll", class_variable)
                print "AIC score with CLL for node " + class_variable + " in the first BN: " + str(aic_2)
                print "AIC increase for node " + class_variable + " :" + str(aic_2 - aic_1)
                
            else:
                #CLL
                cond_1=None
                cond_2=None
                
                #AIC
                aic_1=None
                aic_2=None
                
                
            # Count correct classifications
            count_class1 = self.accuracy(variables, data, class_variable, classification_threshold, bn_1)
            print "Accuracy for the first BN: " + str(count_class1)
            count_class2 = self.accuracy(variables, data, class_variable, classification_threshold, bn_2)
            print "Accuracy for the second BN: " + str(count_class2)
                
        return inverted, added, removed, likelihood_bn1, likelihood_bn2, bic_bn1, bic_bn2, cond_1, cond_2, count_class1, count_class2, aic_1, aic_2
            

    def __parallelConditionalLogLikelihood(self, variables, chunk, class_variable, bn):
        cll = decimal.Decimal('0.0')
        for example in chunk:
            var_val = zip(variables, example)
            # Check if the class variable exist and has a value
            ind_class = variables.index(class_variable)
            if example[ind_class].strip() != "":                    
                class_var_value = var_val.pop(ind_class)
                
                tmp_any = self.any_value
                self.any_value = ""
                query, evidence = self.buildQuery(var_val, [], False)
                self.any_value = tmp_any

                # Try to classify variable
                num = self.variableEliminationInference([class_var_value] + query, [], bn)
                if len(num)==1:
                    num = num[num.keys()[0]]
                elif len(num)==0 :
                    num = decimal.Decimal('0.0')
                
                den = self.variableEliminationInference(query, [], bn)
                if len(den)==1:
                    den = den[den.keys()[0]]
                elif len(den)==0 :
                    den = decimal.Decimal('0.0')
                
                if den != decimal.Decimal('0.0'):
                    inference_result = num/den
                else:
                    inference_result = decimal.Decimal('0.0')
                 
                cll += inference_result.ln() 
                
        return cll       
        
    def conditionalLogLikelihood(self, variables, data, class_variable, bn=None):
        # Computes conditional log likelihood (CLL)
        if bn is None:
            bn=self.BayesianNetwork  
            
        size = data.shape[0]
        n = int(size/self.cores)

        bn = self.getDictRepresentation(bn)
        
        bn_backup = self.BayesianNetwork 
        self.BayesianNetwork = None # workaround for pickle problem

        jobs = []
        job_server.set_ncpus(self.cores)
        if n != 0:
            for i in xrange(0, size, n):
                chunk = data[i:i+n]
                jobs.append(job_server.submit(self.__parallelConditionalLogLikelihood, (variables, chunk, class_variable, bn), modules=('decimal',)))
        else:
            jobs.append(job_server.submit(self.__parallelConditionalLogLikelihood, (variables, chunk, class_variable, bn), modules=('decimal',)))            
        job_server.wait()
        
        total_condlikelihood = 0
        for job in jobs:
            #sys.stdout = open(os.devnull, "w")
            condlikelihood_part = job() # Causes a warning
            #sys.stdout = sys.__stdout__
            total_condlikelihood+=condlikelihood_part
                         
        if self.showTime:
            job_server.print_stats()
        
        self.BayesianNetwork = bn_backup 

        return total_condlikelihood
        
    def expectedConditionalLogLikelihood(self, variables, data, class_variable, ESS, bn=None):
        # Computes conditional log likelihood using ESS to estimate probabilities
        if bn is None:
            bn=self.BayesianNetwork
            
        ecll = decimal.Decimal('0.0')
        for example in data:
            var_val = zip(variables, example)
            # Check if the class variable exist and has a value
            ind_class = variables.index(class_variable)
            if example[ind_class].strip() != "":                    
                class_var_value = var_val.pop(ind_class)
                
                tmp_any = self.any_value
                self.any_value = ""
                query, evidence = self.buildQuery(var_val, [], False)
                self.any_value = tmp_any

                # Try to classify variable
                inference_result = self.inferenceFromProbTable(ESS_all, class_var_value, query)
                if len(inference_result)==1:
                    inference_result = inference_result[inference_result.keys()[0]]
                elif len(inference_result)==0 :
                    inference_result = decimal.Decimal('0.0')
                
                ecll += inference_result.ln()

        return ecll
        
    def mutualInformation(self, node1, node2, bn=None):
        # Computes mutual information for two given nodes.
        if bn is None:
            bn = self.BayesianNetwork
            
        mi = decimal.Decimal('0.0')
        
        for value_n1 in bn[node1].Values:
            for value_n2 in bn[node2].Values:    
                
                # Calculate probability of n1 and n2
                inference_result = self.variableEliminationInference([(node1, value_n1), (node2, value_n2)],[], bn)
                if len(inference_result)==1:
                    inference_result = inference_result[inference_result.keys()[0]]
                elif len(inference_result)==0 :
                    inference_result = decimal.Decimal('0.0')
                    
                prob_n1_n2 = inference_result
                
                # Calculate probability of n1
                inference_result = self.variableEliminationInference([(node1, value_n1)],[], bn)
                if len(inference_result)==1:
                    inference_result = inference_result[inference_result.keys()[0]]
                elif len(inference_result)==0 :
                    inference_result = decimal.Decimal('0.0')
                
                prob_n1 = inference_result
                
                # Calculate probability of n2
                inference_result = self.variableEliminationInference([(node2, value_n2)],[], bn)
                if len(inference_result)==1:
                    inference_result = inference_result[inference_result.keys()[0]]
                elif len(inference_result)==0 :
                    inference_result = decimal.Decimal('0.0')
                
                prob_n2 = inference_result
                
                
                mi += prob_n1_n2 * (prob_n1_n2/(prob_n1)*(prob_n2)).ln()

        return mi
        
    def expectedMutualInformation(self, node1, node2, dataset_size, ESS_all, bn=None):
        # Computes mutual information for two given nodes based on ESS.
        if bn is None:
            bn = self.BayesianNetwork
            
        mi = decimal.Decimal('0.0')
        for value_n1 in bn[node1].Values:
            for value_n2 in bn[node2].Values:    
                # Calculate probability of n1
                cpd_1 = self.eliminateComplementaryLines([(node1, value_n1)], [ESS_all], bn)[0]
                summ_1 = sum(cpd_1.values())
            
                # Calculate probability of n2
                cpd_2 = self.eliminateComplementaryLines([(node2, value_n2)], [ESS_all], bn)[0]
                summ_2 = sum(cpd_2.values())
               
                
                if len(cpd_1) > len(cpd_2):
                    cpd = self.eliminateComplementaryLines([(node1, value_n1), (node2, value_n2)], [cpd_2], bn)[0]
                else:
                    cpd = self.eliminateComplementaryLines([(node1, value_n1), (node2, value_n2)], [cpd_1], bn)[0]
                    
                summ_1_2 =  sum(cpd.values())
                
                if (summ_1*summ_2) != 0: 
                    if summ_1_2 != 0:
                        mi += (summ_1_2/dataset_size) * ((summ_1_2/(summ_1)*(summ_2)) * dataset_size).ln()
        return mi
    
    def BBSL_StepI(self, variables, data, class_variable, classification_threshold, count_threshold, bn, max_nodes):
        bn=self.copy(bn)
        # Structure search - masters degree
        print "Masters search algorithm"
        # Determine the structures to search
        involved_nodes = set([])
        examples_considered = 0
        
        p = ProgressBar()
        count = decimal.Decimal('0.0')
        num_samples = decimal.Decimal(data.shape[0])
        
        # Counter of how many times variables appear in bayes ball result
        bb_count = {}
        for var in variables:
            bb_count[var] = 0
        count_wrong_class = 0
        
        print "     Determining relevant nodes with Bayes Ball..."
        for example in data:
            if self.showProgress:
                p.render(int(count/num_samples * 100), 'Determining relevant nodes.')

            count += decimal.Decimal('1.0')

            inference_result, observed_variables = self.classifyExample(variables, example, class_variable, bn)
            
            if inference_result is not None:
                examples_considered+=1

                # Check if the classification is made wrong
                if inference_result <= classification_threshold:
                    #If so, we need to get the nodes involved using bayes ball
                    observed_variables.remove(class_variable)
                    irrelevant_nodes, requisite_probability_nodes, requisite_observation_nodes, visited_nodes = self.bayesBall([class_variable], observed_variables, bn)
                    relevant_nodes = (set(visited_nodes))
                    # Count how many times each node appears
                    for var in variables:
                        if var in relevant_nodes:
                            bb_count[var] = bb_count[var]+1
                    count_wrong_class += 1
                    
                    involved_nodes = involved_nodes.union(relevant_nodes)
                        
        # Filter nodes that doesn't appear very much in the results
        elim_nodes = []
        for var in involved_nodes:
            perc = float(bb_count[var])/float(count_wrong_class)           
            if perc < count_threshold:
                elim_nodes.append(var)
        # print nodes eliminated
        for var in elim_nodes:
            involved_nodes = involved_nodes.difference(set([var]))
            
            print "Eliminated: " + var + " " + str(float(bb_count[var])/float(count_wrong_class)) + " " + str(bb_count[var]) + "/" + str(count_wrong_class)
            
        # Maximum number of nodes filter:
        if max_nodes is not None:
            sorted_bb_count = sorted(bb_count.iteritems(), key=operator.itemgetter(1))
            sorted_bb_count.reverse()
            to_be_removed = sorted_bb_count[max_nodes:]
            for i in range(len(to_be_removed)):
                node_name = to_be_removed[i][0]
                if node_name in involved_nodes:
                    involved_nodes.remove(node_name)
                    print "Removed " + node_name + " by maximum number of nodes criteria ("+ str(max_nodes) +"). " 
                
        if self.showProgress:
            p.render(100, 'Determining relevant nodes.')
            
        print "     Examples considered: " + str(examples_considered)    
            
        return involved_nodes
            
            
    def BBSL(self, variables, data, class_variable, classification_threshold=0.5, ESS_all=None, score_type='bic', mutual_information_filter=None,  bn=None, external_edges=True, count_threshold=1, involved_nodes=None):
        bn=self.copy(bn)
           
        if involved_nodes is None:  
            involved_nodes = self.BBSL_StepI(variables, data, class_variable, classification_threshold, count_threshold, bn)
          
        # Now that the nodes that can improve classification are determined we may evaluate the possible networks
        # Compute ESS
        if ESS_all is None:
            print "     Calculating all expected sufficient statistics..."
            #ESS_all, combinations = self.computeAllESS(variables, data, bn, None)
            ESS_all, comb1, comb2 = self.computeESS(variables, data, None, bn, False, None, involved_nodes)
            print "     Finished calculating all expected sufficient statistics."

        # Generate possible structures
        print "     " + str(len(involved_nodes)) + " nodes choosen by Bayes Ball."
        
        if external_edges:
            external_edges_to_all = [class_variable]
        else:
            external_edges_to_all = None
            
        best_structure = self.greedyHillClimbing(variables, data, complete_data=False, score_type=score_type, ESS_all=ESS_all, initial_structure=bn, edges_to_consider=None, max_iter=decimal.Decimal("1"), class_variable=class_variable, nodes_to_consider=list(involved_nodes), external_edges_to_all=external_edges_to_all)
            
        print "     Finished evaluating structures."
        
        # Return the best structure found
        return best_structure


    def nFoldCrossValidation(self, num_folds, variables, data, experiment_path, timeout_param=0, partial_out_path=None):
        # Executes an n fold cross validation experiment.
        size = decimal.Decimal(data.shape[0])
        n_per_fold = decimal.Decimal(size/num_folds).to_integral_exact(rounding=decimal.ROUND_FLOOR)
        rest = size%num_folds

        folds =[]
        i = 0
        while i < size:
            in_fold = n_per_fold
            if rest > 0:
                in_fold+=1
                rest-=1
            folds.append(data[(i+1):i+in_fold+1])
            i+=in_fold
           
        # Loading module
        mod = imp.load_source('Experiment', experiment_path)
        exp = mod.Experiment()
                        
        # Leave-one-out
        k = 0                
        while k<num_folds:
            training = copy.deepcopy(folds) 
            test = training.pop(k)
            
            d_training = training[0]
            for w in range(1,len(training)):
                d_training = np.concatenate((d_training,training[w]), axis=0)
            
            # Experiment
            # TRAINING: returns the obtained structure
            #sys.stdout = open(os.devnull, "w")
            @timer.timeout(timeout_param)
            def f():
                training_result_structure = exp.training(variables, d_training)
                return training_result_structure
            
            try:
                training_result_structure = f()
                # TEST: returns the evaluation of the structure
                result = exp.test(variables, test, training_result_structure)
            except timer.TimeoutError:
                if partial_out_path != None:
                    partial_out_path = partial_out_path.replace("{FOLD}", str(k))
                    if os.path.exists(partial_out_path):
                        partial_structure = self.parse(partial_out_path)
                        result = exp.test(variables, test, partial_structure)
                        print "Timeout for this fold. Partial results available."
                    else:
                        result = "Timeout for this fold. Going to the next fold if there's one."
                else:
                    result = "Timeout for this fold. Going to the next fold if there's one."
                
            #sys.stdout = sys.__stdout__
            k+=1
            
            print result
            
    def g2(self, ESS_all, variable1, variable2, conditioned, bn):
        # Conditional Independence Test G2, assimptotically distributed as chi-squared
        # Returns the p-value
        # From "The Max Min Hill Climbing Bayesian Network Structure Learning Algorithm"
        # Under the null hypothesis of conditional independence holding
        # The solution to the zeroes was not described in the article
        
        # Compute G2 statistic value
        g2 = 0
        for comb_cond in self.generateVariableValuesCombinations(conditioned, bn):

            # Check if comb_cond is empty            
            if len(comb_cond) > 0:
                ESS_conditioned = self.eliminateComplementaryLines(comb_cond, [ESS_all], bn)[0]
            else:
                ESS_conditioned = ESS_all    
                
            s_cond = sum(ESS_conditioned.values())
            if s_cond == decimal.Decimal('0'):
                s_cond = self.very_small_value
                
            for value_v1 in bn[variable1].Values:
                for value_v2 in bn[variable2].Values:

                    # Condition on the first variable
                    query = []
                    query.append((variable1, value_v1))
                    ESS_restricted_v1 = self.eliminateComplementaryLines(query, [ESS_conditioned], bn)[0]
                    s_cond_v1 = sum(ESS_restricted_v1.values())
                    if s_cond_v1 == decimal.Decimal('0'):
                        s_cond_v1 = self.very_small_value
                    
                    # Condition on the second variable
                    query = []
                    query.append((variable2, value_v2))
                    ESS_restricted_v2 = self.eliminateComplementaryLines(query, [ESS_conditioned], bn)[0]
                    s_cond_v2 = sum(ESS_restricted_v2.values())
                    if s_cond_v2 == decimal.Decimal('0'):
                        s_cond_v2 = self.very_small_value
                                            
                    # Condition on both variables
                    query.append((variable1, value_v1))
                    ESS_restricted_v1_v2 = self.eliminateComplementaryLines(query, [ESS_conditioned], bn)[0]    
                    s_cond_v1_v2 = sum(ESS_restricted_v1_v2.values())    
                    if s_cond_v1_v2 == decimal.Decimal('0'):
                        s_cond_v1_v2 = self.very_small_value                    
                    
                    # Add to G2 statistic sum
                    g2 += s_cond_v1_v2 * ((s_cond_v1_v2 * s_cond)/(s_cond_v1 * s_cond_v2)).ln()
                    
                    
        g2 = decimal.Decimal("2") * g2      

        # Compute the number of degrees of freedom
        dof = decimal.Decimal((len(bn[variable1].Values) - 1) * (len(bn[variable2].Values) - 1))
        for var_cond in conditioned:
            dof =  dof * decimal.Decimal(len(bn[var_cond].Values))
            
        # Compute p-value
        pvalue =  1 - chi2.cdf(float(g2), float(dof)) 
         
        print "pvalue for G2 statistic for " + variable1 + " and " + variable2 + " conditioned on " + ",".join(conditioned) + " is: " + str(pvalue)     
        return pvalue
    
    def subsets(self, values):
        # Returns all the subsets of a given set of values, expects a list
        subsets = []
        subsets.append([])
        for i in range(len(values)):
            iter = it.combinations(values, i+1)
            for k in iter:
                subsets.append(list(k))
                
        return subsets
    
    def MinAssoc(self, ESS_all, var, target_variable, CPC, bn):
        # Describes the minimum association of var and target_variables given all the subsets of CPC
        # From "The Max Min Hill Climbing Bayesian Network Structure Learning Algorithm"
        all_subsets = self.subsets(CPC)
        min = decimal.Decimal("+inf")
        print "Executing MinAssoc for vars: " + var + " and " + target_variable + ". CPC is : " + ",".join(CPC)
        for a_subset in all_subsets:            
            p_value = self.g2(ESS_all, var, target_variable, a_subset, bn)
            print "p_value for subset: " + ",".join(a_subset) + " is: " + str(p_value)
            
            if p_value < min:
                print "Found minimum value."
                min = p_value
                
        return  min
    
    def MaxMinHeuristic(self, target_variable, CPC, ESS_all, bn):
        # MaxMinHeuristic of Max Min Hill Climbing
        max = decimal.Decimal("-inf")
        print "Executing Max Min Heuristic..."
        for var in bn.keys():
            if var != target_variable:
                min = self.MinAssoc(ESS_all, var, target_variable, CPC, bn)
                min_var = var
                
                if min > max:
                    print "Found maximum value."
                    max = min
                    max_var = min_var
                    
        print "Return from MaxMin Heuristic is value " + str(max) + " from variable " + max_var
        return max_var, max
    
    def MMPC_bar(self, target_variable, ESS_all, bn):
        # Implementation of a step of MMPC 
        # From "The Max Min Hill Climbing Bayesian Network Structure Learning Algorithm" 
        min_p_value = decimal.Decimal("0.05")
        print "Executing MMPC_bar with min p-value of " + str(min_p_value)
        # Phase I: Forward
        print "Entered forward phase."
        CPC = []
        last_CPC = None
        while (last_CPC != CPC):
            last_CPC = copy.deepcopy(CPC)
            
            assoc_var,assoc_value = self.MaxMinHeuristic(target_variable, CPC, ESS_all, bn)
            
            if assoc_value > min_p_value:
                if assoc_var not in CPC:
                    print assoc_var + " added to CPC."
                    CPC.append(assoc_var)
                
        # Phase II: Backward
        print "Entered backward phase"
        all_subsets = self.subsets(CPC)
        final_CPC = []
        for var in CPC:
            for a_subset in all_subsets:
                p_value = self.g2(ESS_all, var, target_variable, a_subset, bn)
                print "p-value for variable " + var + " and target variable " + target_variable + ", conditioned on CPC " + ",".join(a_subset) + " is " + str(p_value)
                if p_value > min_p_value:
                    print "Variable added to MMPC output."
                    final_CPC.append(var)
                    break
        
        return final_CPC
            
    def MMPC(self, target_variable, ESS_all, bn):
        # Implementation of MMPC
        # Determines the minimum set of parents and children of the target_variable
        # From "The Max Min Hill Climbing Bayesian Network Structure Learning Algorithm"
        print "Executing MMPC for variable " + target_variable
        final_CPC = []
        target_CPC = self.MMPC_bar(target_variable, ESS_all, bn)
        print "The target CPC is: " + ",".join(target_CPC)
        for var in target_CPC:
            var_CPC = self.MMPC_bar(var, ESS_all, bn)
            print "The CPC for variable " + var + " is: " + ",".join(var_CPC)
            if target_variable in var_CPC:
                print "Variable included in MMPC result."
                final_CPC.append(var)        
        return final_CPC
    
    def MMHC(self, variables, data, ESS_all, bn):
        # Implementation of MMHC (Max Min Hill Climbing)
        # From "The Max Min Hill Climbing Bayesian Network Structure Learning Algorithm"
        # BN is the initial structure to be considered
        edges_to_consider = []
        print "Executing MMHC..."
        for variable in bn.keys():
            pc = self.MMPC(variable, ESS_all, bn)
            for pc_var in pc:
                edges_to_consider.append((variable, pc_var))       
        
        # Use greedy hill climbing restricted to the set of possible changes already determined in possible_rel
        print "Edges to be considered in GHC: "
        print edges_to_consider
        best_bn = self.greedyHillClimbing(variables, data, False, 'likelihood', ESS_all, bn, edges_to_consider)
        
        return best_bn
                
        
    def markovBlanketFromData(self, ESS_all, query_variable, bn, min_pvalue=decimal.Decimal("0.05")):
        # Identifying the Markov Blanket using the Grow and Shrink algorithm part
        print "Identifying the Markov Blanket for " + query_variable + "..."
        blanket = []
        variables = bn.keys()
        
        # Growing Phase
        variables.remove(query_variable)
        for var in variables:
            # Condtional Independence Test
            pvalue = self.g2(ESS_all, query_variable, var, blanket, bn)
            # If dependent
            if pvalue < min_pvalue:
                blanket.append(var)
        
        # Shrinking Phase
        for var in blanket:
            blanket_without_var = copy.deepcopy(blanket)
            blanket_without_var.remove(var)
            
            #Conditional independence test
            pvalue = self.g2(ESS_all, query_variable, var, blanket_without_var, bn)
            # if independent
            if pvalue > min_pvalue:    
                blanket.remove(var)
                
        return blanket
    
    def computeAllMarkovBlankets(self, ESS_all, bn, min_pvalue=decimal.Decimal("0.05")):
        # Computes the markov blankets of all variables in the bayesian networks from data
        # Step of Grow and Shrink algorithm
        print "Computing Markov Blankets of all nodes..."
        variables = bn.keys()
        
        blankets = {}
        for var in variables:
            blanket = self.markovBlanketFromData(ESS_all, var, bn, min_pvalue)
            blankets[var] = blanket
        
        return blankets
    
    def computeStructureFromBlankets(self, blankets, bn, ESS_all, min_pvalue=decimal.Decimal("0.05")):
        # Build the structure of the graph using the Markov Blankets of the nodes
        print "Building the structure using the blankets..."
        variables = bn.keys()
        bn = self.copy(bn)
        
        direct_neighbours = {}
        for variable_1 in variables:
            direct_neighbours[variable_1] = []
            for variable_2 in blankets[variable_1]:
                dependent = True
                # Choose the set with less elements

                if len(blankets[variable_1]) >= len(blankets[variable_2]):
                    t_set = copy.deepcopy(blankets[variable_1])
                    if variable_2 in t_set:
                        t_set.remove(variable_2)
                else:
                    t_set = copy.deepcopy(blankets[variable_2])
                    if variable_1 in t_set:
                        t_set.remove(variable_1)
                    
                # Generate all subsets    
                all_subsets = self.subsets(t_set)
                for s_subset in all_subsets:
                    # Conditional Independence Tests
                    pvalue = self.g2(ESS_all, variable_1, variable_2, s_subset, bn)
 
                    # if independent
                    if pvalue > min_pvalue:
                        dependent=False
            
                # If there isn't a set that makes X and Y independent, they're direct neighbours         
                if dependent:
                    direct_neighbours[variable_1].append(variable_2)
                    bn[variable_1].Parents.append(variable_2)
                    bn[variable_2].Parents.append(variable_1)

                
        return bn, direct_neighbours
    
    def orientEdgesGS(self, ESS_all, blankets, direct_neighbours, bn, min_pvalue=decimal.Decimal("0.05")):
        # Orient Edges step of Grow and Shrink algortihm 
        # From "Bayesian Network Induction via Local Networks" (1999)
        print "Orienting edges in Grow and Shrink (Step 3)..."
        variables = bn.keys()
        bn= self.copy(bn)
        
        for variable_1 in variables:
            if len(direct_neighbours[variable_1])>0:
                for variable_2 in direct_neighbours[variable_1]:
                    exists_z = False
                    # Check if exists Z
                    z_set = set(direct_neighbours[variable_1]).difference(set(direct_neighbours[variable_2])).difference(set([variable_2]))
                    z_set = list(z_set)
                    
                    for variable_3 in z_set:
                        dependent = True
                        # Check if Z and Y are dependent  given S union {X}
                        # Computing U set
                        b_y = set(blankets[variable_2]).difference(set([variable_1, variable_3]))
                        b_z = set(blankets[variable_3]).difference(set([variable_1, variable_2]))
                        b_y = list(b_y)
                        b_z = list(b_z)
                        # Get the smaller
                        if len(b_y) <= len(b_z):
                            u_set = b_y
                        else:
                            u_set = b_z
                            
                        # Generate subsets
                        all_subsets = self.subsets(u_set)
    
                        for u_subset in all_subsets:
                            # Conditional independence teste
                            # Check if Y and Z are independent given S (u_subset) union X
                            u_subset = u_subset.append(variable_1)
                            pvalue = self.g2(ESS_all, variable_2, variable_3, u_set, bn)
                            
                            if pvalue > min_pvalue:
                                dependent=False
                                
                        if dependent:
                            exists_z = True
                            
                    if exists_z:
                        # Add Y as parent of X 
                        if not (variable_2 in bn[variable_1].Parents):
                            bn[variable_1].Parents.append(variable_2)
                        exists_z = False

        return bn
    
    
    def bnToNetworkx(self, bn):
        # Converts a bayesian network in BN format to NetworkX format
        bn_networkx = networkx.DiGraph()
        # Add nodes
        bn_networkx.add_nodes_from(bn.keys())
        for node in bn.keys():
            parents = bn[node].Parents
            n_parents = len(bn[node].Parents)
            
            # Add edges
            edges = zip(parents, n_parents * [node])
            if len(edges)>0:
                bn_networkx.add_edges_from(edges)
        
        return bn_networkx
    
    def removeCyclesGS(self, bn):
        # Step 4 of Grow and Shrink
        # Remove the cycles as in  "Bayesian Network Induction via Local Networks" (1999)
        print "Removing Cycles in Grow and Shrink (Step 4)..."
        bn = self.copy(bn)
        r_set = []
        cycles = [0] # Just to force entering loop
        while (1):
            # Convert bn to Networkx format
            bn_networkx = self.bnToNetworkx(bn)
                    
            # Compute cycles in the graph
            cycles = networkx.simple_cycles(bn_networkx)
    
            # Build list of edges
            edges = []
            for cycle in cycles:
                # multiply nodes in cycle by its shifted version to get edges in a path
                shift_cycle = copy.deepcopy(cycle)
                shift_cycle.pop(0)
                edges = edges + zip(cycle, shift_cycle)
            
            if len(cycles) > 0:
                # Get edge that participates in most cycles
                count = collections.Counter(edges)
                max_tuple = count.most_common()[0][0]
                
                # Append edge to r_set and remove from graph
                node = max_tuple[1]
                parent = max_tuple[0]
    
                bn[node].Parents.remove(parent)    
                r_set.append(max_tuple)
            else:
                break
        
        return bn, r_set
      
    def reverseEdgesGS(self, bn, r_set):
        # Reverse given edges and insert it into bn
        print "Reversing edges in Grow and Shrink (Step 5)..."
        bn = self.copy(bn)
        for edge in r_set:
            node = edge[1]
            parent = edge[0]
            if not self.existEdge(parent, node, bn):
                # Add inverted edge
                bn[parent].Parents.append(node)
            
        return bn
    
    def propagateDirectionsGS(self, direct_neighbours, bn):
        # Final step of Grow and Shrink algorithm
        # From "Bayesian Network Induction via Local Networks" (1999)
        print "Propagating directions in Grow and Shrink (Step 6)..."
        bn=self.copy(bn)
        for node in bn.keys():
            for neighbour in direct_neighbours[node]:
                # Is there an edge between node and neighbour?
                if not (self.existEdge(node, neighbour, bn)):
                    bn_networkx = self.bnToNetworkx(bn)
                    # Is there a path from node to neighbour now?
                    if networkx.has_path(bn_networkx, node, neighbour):
                        print "Added edge: ( "+node+" , "+neighbour+" )"
                        bn[neighbour].Parents.append(node)
        return bn 
                
    def GrowShrink(self, ESS_all, bn, min_pvalue=decimal.Decimal("0.05")):
        # Implementation of the Grow and Shrink structure learning algorithm
        # From "Bayesian Network Induction via Local Networks" (1999)
        print "Executing Grow and Shrink algorithm..."
        bn=self.copy(bn)
   
        #Compute Markov Blankets
        blankets = self.computeAllMarkovBlankets(ESS_all, bn, min_pvalue)
        
        #ComputeGraphStructure
        bn_not_used, direct_neighbours = self.computeStructureFromBlankets(blankets, bn, ESS_all, min_pvalue)
        #Orient Edges
        bn = self.orientEdgesGS(ESS_all, blankets, direct_neighbours, bn)

        # Remove Cycles
        bn, r_set = self.removeCyclesGS(bn)

        # Reverse Edges
        bn = self.reverseEdgesGS(bn, r_set)
        
        # Propagate Directions
        bn = self.propagateDirectionsGS(direct_neighbours, bn)

        return bn

    def isIndependent(self, node, bn):
        # Checks if a node is independent of all others in the graph
        if len(bn[node].Parents)==0:
            for possible_child in bn.keys():
                if node  in bn[possible_child].Parents:
                    return False
            return True
        else:
            return False
        
    
    def chisquared(self, ESS_all, variable1, variable2, conditioned, bn):       
        # Determines the p-value for the chi-squared statistics with df degrees of freedom
        chi_squared = 0
        for comb_cond in self.generateVariableValuesCombinations(conditioned, bn):
            # Check if comb_cond is empty            
            if len(comb_cond) > 0:
                ESS_conditioned = self.eliminateComplementaryLines(comb_cond, [ESS_all], bn)[0]
            else:
                ESS_conditioned = ESS_all    
                
            s_cond = sum(ESS_conditioned.values())
            if s_cond == decimal.Decimal('0'):
                s_cond = self.very_small_value
                
            for value_v1 in bn[variable1].Values:
                for value_v2 in bn[variable2].Values:

                    # Condition on the first variable
                    query = []
                    query.append((variable1, value_v1))
                    ESS_restricted_v1 = self.eliminateComplementaryLines(query, [ESS_conditioned], bn)[0]
                    s_cond_v1 = sum(ESS_restricted_v1.values())
                    if s_cond_v1 == decimal.Decimal('0'):
                        s_cond_v1 = self.very_small_value
                    
                    # Condition on the second variable
                    query = []
                    query.append((variable2, value_v2))
                    ESS_restricted_v2 = self.eliminateComplementaryLines(query, [ESS_conditioned], bn)[0]
                    s_cond_v2 = sum(ESS_restricted_v2.values())
                    if s_cond_v2 == decimal.Decimal('0'):
                        s_cond_v2 = self.very_small_value
                                            
                    # Condition on both variables
                    query.append((variable1, value_v1))
                    ESS_restricted_v1_v2 = self.eliminateComplementaryLines(query, [ESS_conditioned], bn)[0]    
                    s_cond_v1_v2 = sum(ESS_restricted_v1_v2.values())    
                    if s_cond_v1_v2 == decimal.Decimal('0'):
                        s_cond_v1_v2 = self.very_small_value                    
                    
                    # Add to chi-squared statistic sum
                    chi_squared += (((s_cond * s_cond_v1_v2) - (s_cond_v1 * s_cond_v2))**2)/(s_cond * s_cond_v1 * s_cond_v2)
                    

        # Compute the number of degrees of freedom
        dof = decimal.Decimal((len(bn[variable1].Values) - 1) * (len(bn[variable2].Values) - 1))
        for var_cond in conditioned:
            dof =  dof * decimal.Decimal(len(bn[var_cond].Values))
            
        # Compute p-value
        pvalue =  1 - chi2.cdf(float(chi_squared), float(dof)) 
         
        print "pvalue for Chi-Squared statistic for " + variable1 + " and " + variable2 + " conditioned on " + ",".join(conditioned) + " is: " + str(pvalue)     
        return pvalue 
        
    def k2_score(self, ESS_all, node, parents, bn):
        #Computes the log of the score used in k2 structure search algorithm
        # From: Illustration of the K2 Algorithm for Learning Bayes Net Structures, Carolina Ruiz,Department of Computer Science, WPI
        vars = bn.keys()
        
        # Calculate log of one of the terms
        alpha = decimal.Decimal("0.0")
        for comb in self.generateVariableValuesCombinations(parents + [node], bn):
            ESS_conditioned = self.eliminateComplementaryLines(comb, [ESS_all], bn)[0]
            summ = sum(ESS_conditioned.values())
            fat = math.factorial(round(summ))
            alpha = decimal.Decimal(str(fat)).ln() + alpha
                    
            
        # Calculate log of the other term
        num_values = len(bn[node].Values)
        for comb in self.generateVariableValuesCombinations(parents, bn):
            numerator = decimal.Decimal(math.factorial(num_values - 1))
            ESS_conditioned = self.eliminateComplementaryLines(comb, [ESS_all], bn)[0]
            nij = sum(ESS_conditioned.values()) 
            denominator = decimal.Decimal(math.factorial(round(nij) + num_values - 1))
            first_term = (numerator/denominator).ln()
        
        return first_term * alpha
        
    
    def k2(self, variables, data, node_order, ESS_all, max_parents, bn):
        # K2 algortithm for structure learning
        parents_node = {}
        node_order.reverse()
        
        if ESS_all is None:
            ESS_all = self.computeAllESS(variables, data, bn, None, False, [])
        
        predecessors = []
        for node in node_order:
            print "Node: " + node
            parents_node[node] = []
            # Calculate score the first time
            old_score = self.k2_score(ESS_all, node, [], bn)
            print "Old Score: " + str(old_score) 
            proceed = True
            while (proceed) and (len(parents_node[node]) < max_parents):
                new_score = decimal.Decimal("-inf")
                predecessors_without_parents = list(set(predecessors).difference(set(parents_node[node])))
                print "Predecessors (without parents): " + ",".join(predecessors_without_parents)
                for predecessor in predecessors_without_parents:
                    tmp_score = self.k2_score(ESS_all, node, parents_node[node] + [predecessor], bn)
                    print "Predecessor: " + predecessor + ", score: " + str(tmp_score)
                    if  (tmp_score > new_score):
                        new_score = tmp_score
                        z = predecessor
                       
                if new_score > old_score:
                    old_score = new_score
                    parents_node[node].append(z)
                    print "Node: " + node + " Parent: " + z
                else:
                    proceed = False
                    
            predecessors.append(node)
            
        for node in parents_node:
            bn[node].Parents = parents_node[node]  
            
        size = data.shape[0]                      
        bn = self.fixedParameterBNInitialization(decimal.Decimal('0.0'), bn, [], True)
        bn = self.estimateBNParamsFromESS(ESS_all, bn, bn.keys())       

        return bn
    
    def randomPerturbation(self,structure, ESS_all, node=None,  distance=None, use_markov_blanket=False):
        # Generates a network by making a random change to the actual network
        #Define which change is gonna be made
        
        # IF LOCAL
        # Generates a network by making a local change to the given network.
        # Changes are made to the adjacency of the node. 
        # Changes are random and made to the edges in the markov blanket, or made to the edges that join nodes in a given distance from the main node.
        
        if node is not None:
            # Determine edges to consider
            # LOCAL APPROACH
            if use_markov_blanket==True:
                nodes_to_choose = self.getMarkovBlanket(node, structure)
            
            if distance is not None:
                distance_bfs, nodes_to_choose = self.BFS(structure, node, None, distance)
        else:
            # GLOBAL APPROACH
            nodes_to_choose = structure.keys()
        
        if len(nodes_to_choose) >=2:
            # Choose two nodes randomly
            change_achieved = False
            while not change_achieved: 
                first_var = random.choice(nodes_to_choose)
                nodes_to_choose.remove(first_var)
                second_var = random.choice(nodes_to_choose)
                changed_connection = [first_var, second_var]
                
                # Check connection status
                exist_edge = self.existEdge(first_var, second_var, structure)
                if exist_edge:
                    operators = ["remove", "invert"]
                else:
                    operators=["add", "invert"]
                
                # Choose operator
                apply_operator = random.choice(operators)
                operators.remove(apply_operator)
                other_operator = operators[0]
                
                # Now, apply the operator to the network and check for aciclicity
                bn = self.copy(structure)
                if apply_operator=="invert":
                    if not exist_edge:
                        bn = self.addEdge(first_var, second_var, bn)
                    bn = self.invertEdge(first_var, second_var, bn)
                if apply_operator=="add":
                    bn = self.addEdge(first_var, second_var, bn)
                if apply_operator=="remove":
                    bn = self.removeEdge(first_var, second_var, bn)
                 
                is_aciclic = self.checkAciclicity(bn)
                if not is_aciclic:
                    #Apply other operator
                    bn = self.copy(structure)
                    if other_operator=="invert":
                        if not exist_edge:
                            bn = self.addEdge(first_var, second_var, bn)
                        bn = self.invertEdge(first_var, second_var, bn)
                    if other_operator=="add":
                        bn = self.addEdge(first_var, second_var, bn)
                    if other_operator=="remove":
                        bn = self.removeEdge(first_var, second_var, bn)            
                    is_aciclic = self.checkAciclicity(bn)
                    if not is_aciclic:
                        bn = self.copy(structure)
                        #print "RandomPerturbation: Network not changed. Trying again."
                    else:
                        change_achieved = True
                else:
                    change_achieved = True
                    
            bn = self.fixedParameterBNInitialization(decimal.Decimal("0.0"), bn, changed_connection, True)          
            bn = self.estimateBNParamsFromESS(ESS_all, bn,  changed_connection)
        else:
            print "Not sufficient nodes. Less than two nodes, need at least 2 to complete."
            bn = None    
        return bn
    
    def adjacencies(self, bn, search_node=None):
        # Returns the adjacencies of a node.
        # Also loads all adjacencies of all nodes in self.adjacencies_dict
        
        if self.adjacencies_dict is None:
            #Initialize dictionary
            adjacencies_dict = {}
            for node in bn.keys():
                adjacencies_dict[node] = [] 
            
            #Build adjacencies dictionary
            for node in bn.keys():
                for parent in bn[node].Parents:
                    adjacencies_dict[node].append(parent)
                    adjacencies_dict[parent].append(node)
            
            self.adjacencies_dict = adjacencies_dict
                
            if search_node is None:
                return adjacencies_dict
            else:
                return adjacencies_dict[search_node]
        else:
            return self.adjacencies_dict[search_node]
    
    
    def BFS(self, structure, start_node, search_node, depth_limit):
        # Runs a breadth first search and returns the distance to the search node, and also the visited nodes.
        queue = []
        queue_distance = []
        marked = []
        queue.append(start_node)
        queue_distance.append(0)
        marked.append(start_node)
        
        while (len(queue) != 0):
            node = queue.pop(0)
            distance = queue_distance.pop(0)
            
            if node == search_node:
                return distance, marked
            
            if distance < depth_limit:

                adjacencies = self.adjacencies(structure, node)
                for adj_node in adjacencies:
                    if adj_node not in marked:
                        marked.append(adj_node)
                        queue.append(adj_node)
                        queue_distance.append(distance+1)
        
        distance = -1
        return distance, marked      
       
    
    def simulatedAnnealing(self, variables, data, ESS_all, ini_structure, perturbation_func, class_variable=None, ini_temp=decimal.Decimal("1000"), temp_threshold=decimal.Decimal("10"), decreasing_temp_factor=decimal.Decimal("0.70"),  score_type='likelihood', max_iter=3):
        # Simulated Annealing algorithm implementation
        # Extracted from Parallel Two-Level Simulated Annealing
        print "Started Simulated Annealing..."
        if ini_structure is None:
            ini_structure = self.BayesianNetwork
            
        if score_type=='likelihood':
            f_old = self.logLikelihood(variables, data, ini_structure)
        elif score_type=='cll':
            f_old = self.conditionalLogLikelihood(variables, data, class_variable, ini_structure)
        
        x_old = self.copy(ini_structure)
        x_best = self.copy(ini_structure)
        f_best = f_old
        print "Initial score:" + str(f_old)
        temperature = ini_temp
        print "Temperature: " + str(temperature)
        while (temperature > temp_threshold):
            for i in range(max_iter):
                st = time.time()
                x_new = perturbation_func(x_old, ESS_all)
                if score_type=='likelihood':
                    f_new = self.logLikelihood(variables, data, x_new)
                elif score_type=='cll':
                    f_new = self.conditionalLogLikelihood(variables, data, class_variable, x_new)
                r = random.random()
                
                print ""
                print "Score: " + str(f_new) + " Time: " + str(time.time()-st)
                boltzman = math.exp((f_new - f_old)/temperature)
                print "Boltzman constant: " + str(boltzman)
                
                if (f_new>=f_old) or (r<=boltzman):
                    if (f_new<f_old) and (r<=boltzman):
                        print "Worse but changed network."
                    x_old = self.copy(x_new)
                    f_old = f_new
                    if f_new > f_best:
                        x_best = self.copy(x_new)
                        f_best = f_new
                        print "Best score: " + str(f_best)
                        
            temperature = decreasing_temp_factor * temperature
            print ""
            print "New temperature: " + str(temperature)
        
        print "SA Best Score: " + str(f_best)    
        return x_best, f_best    
        
        
    def TwoLevelSimulatedAnnealing(self, variables, data, ESS_all, ini_structure, perturbation_func, optimization_func, class_variable=None, ini_temp=decimal.Decimal("1000"), temp_threshold=decimal.Decimal("10"), decreasing_temp_factor=decimal.Decimal("0.70"),  score_type='likelihood', max_iter=3):
        # Two Level Simulated Annealing algorithm implementation
        # Extracted from Parallel Two-Level Simulated Annealing
        print "Started Two Level Simulated Annealing..."
        if ini_structure is None:
            ini_structure = self.BayesianNetwork
        
        x_old = optimization_func(self.copy(ini_structure), ESS_all)
            
        if score_type=='likelihood':
            f_old = self.logLikelihood(variables, data, ini_structure)
        elif score_type=='cll':
            f_old = self.conditionalLogLikelihood(variables, data, class_variable, ini_structure)
        
        x_best = self.copy(ini_structure)
        f_best = f_old
        print "Initial score:" + str(f_old)
        temperature = ini_temp
        print "Temperature: " + str(temperature)
        while (temperature > temp_threshold):
            for i in range(max_iter):
                st = time.time()
                x_new = perturbation_func(x_old, ESS_all)
                x_new = optimization_func(self.copy(x_new), ESS_all)
                
                if score_type=='likelihood':
                    f_new = self.logLikelihood(variables, data, x_new)
                elif score_type=='cll':
                    f_new = self.conditionalLogLikelihood(variables, data, class_variable, x_new)
                    
                r = random.random()
                
                print ""
                print "Score: " + str(f_new) + " Time: " + str(time.time()-st)
                boltzman = math.exp((f_new - f_old)/temperature)
                print "Boltzman constant: " + str(boltzman)
                
                if (f_new>=f_old) or (r<=boltzman):
                    if (f_new<f_old) and (r<=boltzman):
                        print "Worse but changed network."
                    x_old = self.copy(x_new)
                    f_old = f_new
                    if f_new > f_best:
                        x_best = self.copy(x_new)
                        f_best = f_new
                        print "Best score: " + str(f_best)
                        
            temperature = decreasing_temp_factor * temperature
            print ""
            print "New temperature: " + str(temperature)
        
        print "TLSA Best Score: " + str(f_best)    
        return x_best, f_best   
    
    def classifyExample(self, variables, example, class_variable, bn):
        # Get nodes that are observed and unobserved
        observed_variables = []
        observed_variables_values = []
        unobserved_variables = []
        unobserved_variables_values = []
        for k in range(len(variables)):
            if example[k].strip() != "":
                observed_variables.append(variables[k])
                observed_variables_values.append((variables[k], example[k]))
            else:
                unobserved_variables.append(variables[k])
                unobserved_variables_values.append((variables[k], example[k]))

        # Check if class_variable is observed
        if class_variable in unobserved_variables:
            return None, None
          
        # Get values of the nodes that are not the class_variable
        nodes_values = []
        for tpl in observed_variables_values:
            if tpl[0] == class_variable:
                class_var_value = tpl
            else:
                nodes_values.append(tpl)
                
        # Try to classify variable
        num = self.variableEliminationInference([class_var_value] + nodes_values, [], bn)
        if len(num)==1:
            num = num[num.keys()[0]]
        elif len(num)==0 :
            num = decimal.Decimal('0.0')
        
        if len(nodes_values)!=0:
            den = self.variableEliminationInference(nodes_values, [], bn)
            if len(den)==1:
                den = den[den.keys()[0]]
            elif len(den)==0 :
                den = decimal.Decimal('0.0')
         
            if den != decimal.Decimal('0.0'):
                inference_result = num/den
            else:
                inference_result = decimal.Decimal('0.0')
        else:
            inference_result = num
            
        return inference_result, observed_variables
    
    def countCorrectClassifications(self, variables, data, class_variable, classification_threshold, bn):
        # Count how many classifications are done correctly.
        bn= self.copy(bn)

        examples_considered = 0
        count_right_class = 0
        
        for example in data:
            inference_result, observed_variables = self.classifyExample(variables, example, class_variable, bn)
            
            if inference_result is not None:
                examples_considered+=1
                # Check if the classification is made wrong
                if inference_result > classification_threshold:     
                    count_right_class += 1

        return count_right_class, examples_considered      
 
    def accuracy(self, variables, data, class_variable, classification_threshold=decimal.Decimal("0.5"), bn=None):
        
        size = decimal.Decimal(data.shape[0])
        correct, examples = self.countCorrectClassifications(variables, data, class_variable, classification_threshold, bn)
        
        accuracy = correct/size
        
        return accuracy
        
    def laplaceSmoothBNParams(self, bn, alfa):
        bn = self.copy(bn)
        # Applies laplace smoothing to all the parameters of the bayesian network
        for node in bn.keys():
            CPD = bn[node].CPD
            d = len(bn[node].Values)
            for line in CPD:
                new_parameter = (CPD[line] + alfa)/(decimal.Decimal("1.0") + alfa * decimal.Decimal(d))
                bn[node].CPD[line] = new_parameter

        return bn

    

    def AIC(self, variables, data, bn, score_type="likelihood", class_variable=None):
        # Computes Akaike Information Criterion for a bayesian struture and data
        # score_type can be likelihood or conditional log likelihood
        
        if score_type=="likelihood":
            score = self.logLikelihood(variables, data, bn)
        elif score_type=="cll":
            score = self.conditionalLogLikelihood(variables, data, class_variable, bn)
            
        num_parameters = decimal.Decimal(self.numberOfParameters(bn))
        
        aic = (2*num_parameters) - (2 * score)
        
        return aic

    def pairedTtest(self, values_a, values_b):
        # Returns the pvalue for a paired ttest
        if len(values_a) != len(values_b):
            print "Number of values in the two given sets is different."
            return None
        
        pvalue = ttest_rel(values_a, values_b, 0)[1]
        
        return pvalue
    
    def MB_recursive(self, class_variable, observed_variables, bn):
        # Returns the MB* of the node.
        stack = [class_variable]
        visited = []
        MB_recursive = set([class_variable])
        
        while (len(stack)!=0):
            node = stack.pop(0)
            visited.append(node)
            
            MB_missing = self.getMarkovBlanket(node, bn)
            MB_recursive = MB_recursive.union(set(MB_missing))

            for var in MB_missing:
                if var not in observed_variables:
                    if (var not in visited) and (var not in stack):
                        stack.append(var)

        return set(MB_recursive)
        
    
    def DAHVI(self, variables, data, score, class_variable, classification_threshold, bn):
        # DAHVI algorithm implementation
        # Extracted from  Revisando Redes Bayesianas atraves da Introdução de Variáveis não-observadas
        # By now, it works for only one variable
        
        # Receives initial networ and calculates the score

        current_network = self.copy(bn)
        if (score=='likelihood'):
            current_score = self.logLikelihood(variables, data, current_network)

        if (score=='cll'):
            current_score = self.conditionalLogLikelihood(variables, data, class_variable, current_network)
                              
        if (score == "aic"):
            current_score = self.AIC(variables, data, current_network, "cll", class_variable)  
           
        best_network = current_network 
        best_score=current_score
        best_data = np.copy(data)

        bn_backup = self.BayesianNetwork
        self.BayesianNetwork = None   
        
        progress = True 
        while(progress):
            current_network = best_network
            current_data = best_data
            
            # Catch examples where at least one class variable was wrongly classified (threshold?)
            rev_points=set([])

            for example in current_data:            
                inference_result, observed_variables = self.classifyExample(variables, example, class_variable, current_network)
                
                if inference_result is not None:  
                    # Check if the classification is made wrong
                    if inference_result <= classification_threshold:
                        # Get the MB* of the set of nodes in the structure
                        # Select revision points by union of MB* of each example
                        if (len(observed_variables)<=len(variables)):
                                MB_recursive = self.MB_recursive(class_variable, observed_variables, bn)
                                rev_points = rev_points.union(MB_recursive)
                                
            # if rev_points is empty, nothing to do in the network            
            if len(rev_points)==0:
                return best_network, best_score, best_data
            
            print "Revision Points: " + ",".join(list(rev_points))
            # Generate net neighbours and evaluates it
            progress = False
            for neighbour_network in self.generateDAHVIneighbours(rev_points, ['yes', 'no'], current_network):
                current_data = np.copy(best_data)

                #Add column to the data
                new_vars = list((set(neighbour_network.keys())).difference(set(bn.keys())))
                num_new_cols = len(new_vars)
                size = current_data.shape[0]
                new_cols = np.array([['']*num_new_cols]*size)
                current_data = np.hstack((current_data, new_cols))
                current_variables = variables + new_vars
                
                
                ESS, neighbour_network = self.expectationMaximization(current_variables, current_data, decimal.Decimal('0.01'), neighbour_network, param_ini=None, max_iter=1, local=[])
                del ESS
                print ""
                
                #Evaluate network
                if (score=='likelihood'):
                    current_score = self.logLikelihood(current_variables, current_data, neighbour_network)
                    criteria = (current_score > best_score)
                    
                if (score=='cll'):
                    current_score = self.conditionalLogLikelihood(current_variables, current_data, class_variable, neighbour_network)                  
                    criteria = (current_score > best_score)
                    
                if (score == "aic"):
                    current_score = self.AIC(current_variables, current_data, neighbour_network, "cll", class_variable)    
                    criteria = (current_score < best_score)
                    
                if (criteria):
                    best_network = self.copy(neighbour_network)
                    best_score = current_score
                    best_data = np.copy(current_data)
                    progress = True
                    
                

        self.BayesianNetwork = bn_backup
        return best_network, best_score, best_data
        
        
        
        
        
        
        