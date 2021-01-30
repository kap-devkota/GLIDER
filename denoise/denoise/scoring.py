import random

def score_cv(test_nodes, test_labelling, real_labelling, debug_params = {"DEBUG" : False}):
    """Scores cross validation by counting the number of test nodes that
    were accurately labeled after their removal from the true
    labelling.
    """
    correct = 0
    total   = 0
    if debug_params["DEBUG"] == True:
        if debug_params["D_TYPE"] == "degree_check":
            degree_vec   = debug_params["degree_vec"]
            degree_thres = debug_params["degree_check"]["threshold"]
            n_below      = 0
            n_above      = 0
            c_below      = 0
            c_above      = 0
    for node in test_nodes:
        if node not in test_labelling:
            continue
        if debug_params["DEBUG"] == True:
            if debug_params["D_TYPE"] == "degree_check":
                if degree_vec[node] >= degree_thres:
                    n_above += 1
                else:
                    n_below += 1
        test_label = test_labelling[node]
        if type(test_label) is list:
            for tl in test_label:
                if tl in real_labelling[node]:
                    if debug_params["DEBUG"] == True:
                        if debug_params["D_TYPE"] == "degree_check":
                            if degree_vec[node] >= degree_thres:
                                c_above += 1
                            else:
                                c_below += 1
                    correct += 1
                    break
        else:
            if test_label in real_labelling[node]:
                if debug_params["DEBUG"] == True:
                    if debug_params["D_TYPE"] == "degree_check":
                        if degree_vec[node] >= degree_thres:
                            c_above += 1
                        else:
                            c_below += 1
                correct += 1
        total += 1
    if debug_params["DEBUG"] == True:
        if debug_params["D_TYPE"] == "degree_check":
            b_acc = 0 if n_below == 0 else c_below / n_below
            a_acc = 0 if n_above == 0 else c_above / n_above
            print(f"[!] Error rate for degree below {degree_thres}: {c_below}, {n_below}, {b_acc}")
            print(f"[!] Error rate for degree above {degree_thres}: {c_above}, {n_above}, {a_acc}")
    return float(correct) / float(total)

def kfoldcv(k, labels, prediction_algorithm, randomized=True, degree_vec = None, reverse = False):
    """Performs k-fold cross validation.

    Args:
      - A number of folds k
      - A labeling for the nodes.
      - An algorithm that takes the training labels
      and outputs a predicted labelling.

    Output: 
      - A list where each element is the accuracy of the
      learning algorithm holding out one fold of the data.
    """
    debug_params = {"DEBUG"        : True if degree_vec is not None else False,
                    "D_TYPE"       : "degree_check",
                    "degree_vec"   : degree_vec,
                    "degree_check" : {"threshold" : 10}}
    nodes = list(labels.keys())
    if randomized:
        random.shuffle(nodes)
    accuracies = []
    for i in range(0, k):
        inc = int(len(nodes) / k)
        x = inc * i
        y = inc * (i + 1)
        if i + 1 == k:
            y = len(nodes)
        if not reverse:
            training_nodes = nodes[:x] + nodes[y:]
            test_nodes = nodes[x:y]
        else:
            training_nodes  = nodes[x:y]
            test_nodes      = nodes[:x] + nodes[y:]
            
        training_labels = {n: labels[n] for n in training_nodes}
        test_labelling  = prediction_algorithm(training_labels)
        
        accuracy = score_cv(test_nodes, 
                            test_labelling, 
                            labels, 
                            debug_params = debug_params)
        accuracies.append(accuracy)
    return accuracies
