import itertools
import numpy as np
from collections import defaultdict
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression


def vote(voters, labels_f, weight_f):
    """
    Votes for the most popular label among the voters,
    weighted by their significance. 
    
    Input:
      - A list of voters.
      - A function mapping each voter to a list of labels.
      - A function mapping each voter to its weight.
    Output:
      - The most popoular label or none if no voters have labels.
    """
    label_counts = defaultdict(int)
    for voter in voters:
        for label in labels_f(voter):
            label_counts[label] += weight_f(voter)

    if not label_counts:
        return None

    return max(label_counts.keys(), key=lambda k: label_counts[k])

def wmv(A, labels_f, weight_f=lambda x: x, default_label="????", params = None):
    """
    Weighted majority vote algorithm for an undirected graph.

    Input:
      - An adjacency matrix for a graph.
      - A function mapping node IDs to a list of labels. An
      empty list represents no known label.
      - A function mapping weights to new values. 
      - A label to give when no label is predicted
    Output:
      - A dictionary mapping node IDs to a label. If the label
      is already known, the first label in the list is picked.
    """
    predicted_labels = {}
    
    filter_type = None
    if params  != None:
        filter_type = params["filter_type"]
        if filter_type == "degree":
            deg_vec    = params["degree_vec"]
            max_degree = params["degree_threshold"]            

    n = A.shape[0]
    for i in range(n):
        labels = labels_f(i)
        if labels:
            predicted_labels[i] = labels[0]
            continue
        voters_ = filter(lambda j: A[i, j] != 0, itertools.chain(range(0, i), range(i + 1, n)))
        if filter_type == "degree":
            voters     = []
            count      = 0
            for v in voters_:
                if deg_vec[v] <= max_degree:
                    voters.append(v)
        else:
            voters = voters_
        prediction = vote(voters, labels_f, lambda voter: weight_f(A[i, voter]))
        if prediction is not None:
            predicted_labels[i] = prediction
        else:
            predicted_labels[i] = default_label

    return predicted_labels

def mv(A, labels_f, default_label="????"):
    """
    Unweighted majority vote algorithm for an undirected graph.

    Input:
      - An adjacency matrix for a graph.
      - A function mapping node IDs to a list of labels. An
      empty list represents no known label.
      - A label to give when no label is predicted
    Output:
      - A dictionary mapping node IDs to a label. If the label
      is already known, the first label in the list is picked.
    """
    return wmv(A, labels_f, weight_f=lambda _: 1, default_label=default_label)


def glide(node_association, labels_f, params = None):
    n = len(node_association)
    predicted_labels = {}
    for i in range(n):
        labels = labels_f(i)
        if labels:
            predicted_labels[i] = labels[0]
            continue
        voters = [p for p in node_association[i]]
        weight_f = lambda voter: node_association[i][voter]
        prediction = vote(voters, labels_f, weight_f)
        if prediction is not None:
            predicted_labels[i] = prediction
        else:
            predicted_labels[i] = "????"
    return predicted_labels

def knn(distances, labels_f, k, default_label="????", is_weighted=True, params = None):
    """Performs k-nearest neighors voting algorithm using the passed in
    distance matrix.

    Input:
      - An n x n matrix where each entry represents the distance
      between two points.
      - A function mapping node IDs to a list of labels. An
      empty list represents no known label.
      - A label to give when no label is predicted
    Output:
      - A dictionary mapping node IDs to a label. If the label
      is already known, the first label in the list is picked."""
    predicted_labels = {}
    filter_type      = None

    if params != None:
        filter_type = params["filter_type"]
        if filter_type == "degree":
            deg_vec    = params["degree_vec"]
            max_degree = params["degree_threshold"]
            
        
    n = distances.shape[0] # number of nodes
    for i in range(n):
        labels = labels_f(i)
        if labels:
            predicted_labels[i] = labels[0]
            continue
        voters_ = np.argsort(distances[i, :])
        if filter_type == "degree":
            voters     = []
            count      = 0
            for v in voters_:
                if deg_vec[v] <= max_degree:
                    voters.append(v)
                    count += 1
                if count >= k:
                    break
        else:
            voters = voters_[1:k+1]

        if is_weighted == True:
            weight_f = lambda voter: 1 / (distances[voter, i] if distances[voter, i] > 0 else 0.000001)
        else:
            weight_f = lambda voter: 1            
        prediction = vote(voters, labels_f, weight_f)
        if prediction is not None:
            predicted_labels[i] = prediction
        else:
            predicted_labels[i] = default_label
    return predicted_labels
           
def svm(embedding, labels_f, inv_labels_f = lambda x: x, default_label = "????"):
    """
    Performs SVM classification
    labels_f: returns a class label
    """    
    clf = SVC(gamma = "auto")
    n = embedding.shape[0]
    training = []
    testing  = []
    labels   = []
    for i in range(n):
        l        = labels_f(i) 
        if l:
            training.append(i)
        else:
            testing.append(i)
        labels.append(l)
    clf.fit(embedding[training], np.array(labels[training]))
    predictions = clf.predict(embedding[testing])
    for i, t in enumerate(testing):
        labels[t] = predictions[i]
    return {i : inv_labels_f(j) for i, j in enumerate(labels)}


def jaccard_filter(labels_dct, threshold=0.1):
    """
    filters the set of labels so that no two labels have a jacard similarity
    greater than the threshold
    @param labels_dct: A dictionary with labels as keys and a list of indices
        as its values
    @param threshold: the maximum tolerable similarity
    @return (used_label, unused_labels): the dict passed in split into used and
        unused
    """
    
    keys = list(labels_dct.keys())
    used_labels = {}
    unused_labels = {}

    i = 0
    while i  < len(keys):
        label1 = keys[i]
        proteins1 = labels_dct[label1]
        pro1 = set(proteins1) # used for set intersection and union below
        
        # if we see it here then there is nothing similar to it that we've seen
        # before, and the inner loop removes the similar things coming after
        used_labels[label1] = proteins1
        
        # update the indices 
        j = i + 1
        while j < len(keys):
            label2 = keys[j]
            proteins2 = labels_dct[label2]
            pro2 = set(proteins2)
            
            # compute the jacard index as |inter| / |union|
            jaccard = len(pro1.intersection(pro2)) / len(pro1.union(pro2))
            if jaccard > threshold:
                # too similar to label1 so we cull it
                unused_labels[label2] = proteins2

                # trim out this index for speed-up (note don't increment j)
                keys.pop(j)
            else:
                # using while so we need to manually increment
                j += 1

        # using while so we need to manually increment
        i += 1
    
    return (used_labels, unused_labels)


def jaccard_filter_added_unused(labels_dct, threshold=0.1):
    """
    filters the set of labels so that no two labels have a jacard similarity
    greater than the threshold
    @param labels_dct: A dictionary with labels as keys and a list of indices
        as its values
    @param threshold: the maximum tolerable similarity
    @return (used_label, unused_labels): the dict passed in split into used and
        unused
    """
    
    keys = list(labels_dct.keys())
    used_labels = {}
    unused_labels = {}

    i = 0
    associated_labels = {}
    while i  < len(keys):
        label1 = keys[i]
        associated_labels[label1] = []
        proteins1 = labels_dct[label1]
        pro1 = set(proteins1) # used for set intersection and union below
        
        # if we see it here then there is nothing similar to it that we've seen
        # before, and the inner loop removes the similar things coming after
        used_labels[label1] = proteins1
        
        # update the indices 
        j = i + 1
        while j < len(keys):
            label2 = keys[j]
            proteins2 = labels_dct[label2]
            pro2 = set(proteins2)
            
            # compute the jacard index as |inter| / |union|
            jaccard = len(pro1.intersection(pro2)) / len(pro1.union(pro2))
            if jaccard > threshold:
                # too similar to label1 so we cull it
                associated_labels[label1].append(label2)
                used_labels[label1] += proteins2
                # trim out this index for speed-up (note don't increment j)
                keys.pop(j)
            else:
                # using while so we need to manually increment
                j += 1

        # using while so we need to manually increment
        i += 1
    
    return (used_labels, associated_labels)



def perform_binary_OVA(E, labels, params = {}, clf_type="LR"):
    """
    Perform binary svc on embedding and return the new labels
    @param E: Embedding of size n x k
    @param labels: A dictionary that maps the index in the row of embedding to labels. An index can have many labels
    @param params:
    @return labels: Since the dictionary labels is incomplete (some of the indices donot have any labels associated with it), this function performs SVC for each labels and completels the labels dictionary, and returns it.
    """
    def convert_labels_to_dict(lls):
        """
        This function takes in a list of labels associated with a protein embedding, and returns the dictionary that is keyed .
        by the index of protein embedding with the value 
        """
        l_dct = {}
        for k in lls:
            ll       = lls[k]
            l_dct[k] = {i: True for i in ll}
        return l_dct

    labels_dct = convert_labels_to_dict(labels)
    
    def transpose_labels(labels):
        """
        Returns a dict with go labels as keys with values being a list of 
        proteins with that label
        """
        transpose = {}
        for protein in labels:
            lls = labels[protein]
            for ll in lls:
                if not ll in transpose:
                    transpose[ll] = []
                transpose[ll].append(protein)
        return transpose

    # perform a filter on the label classes we are going to use
    t_labels                     = transpose_labels(labels)
    print(f"The number of All the Labels {len(t_labels)}")

    """
    #[used_labels, unused_labels] = jaccard_filter(t_labels, 0.1)
    [used_labels, assoc_dict]   = jaccard_filter_added_unused(t_labels, threshold=thres)
    print(f"The number of Used Labels {len(used_labels)}")
    """

    samples    = {}
    n          = E.shape[0]

    # Adding Positive samples
    for i in labels:
        lls = labels[i]
        for ll in lls:
            """
            # ignore the labels that we aren't considering
            if ll not in used_labels:
                continue
            """
            if ll not in samples:
                samples[ll] = {"positive" : [], "negative" : [], "null" : [], "clf": None}
                if clf_type != "LR":
                    samples[ll]["clf"] = SVC(gamma = "auto", probability=True)
                else:
                    samples[ll]["clf"] = LogisticRegression(random_state = 0)
            samples[ll]["positive"].append(i)

    # Adding Negative samples and creating null set (unlabeled data)
    null_set = []
    for i in range(n):
        if i not in labels:
            null_set.append(i)
        else:
            for j in samples:
                if j not in labels_dct[i]:
                    samples[j]["negative"].append(i)
    null_set = np.array(null_set)

    # Balance negative and positive samples
    # and train the probabilistic SVMs
    for s in samples:
        n_pos = len(samples[s]["positive"])
        n_neg = len(samples[s]["negative"])
        n_val = n_pos if n_pos < n_neg else n_neg
        samples[s]["positive"] = np.array(samples[s]["positive"][:n_val])
        samples[s]["negative"] = np.array(samples[s]["negative"][:n_val])
        lbls                   = np.zeros((2 * n_val, ))
        lbls[:n_pos]           = 1
        inputs                 = np.concatenate([samples[s]["positive"],
                                                 samples[s]["negative"]])
        samples[s]["clf"].fit(E[inputs], lbls)
    
    print("Here")
    # iterate over the classes computing the probability for each point in
    # the null set for each class
    probabilities = np.zeros(( len(null_set), len(samples) ))
    sample_keys = list(samples.keys())
    for i, s in enumerate(sample_keys):
        # predict_proba returns a matrix, each row is the prediction for a datapoint
        # and each column is for a different class. col 0 is negative, col 1 is positive
        probabilities[:,i] = samples[s]["clf"].predict_proba(E[null_set])[:,1]

    # predict the highest probability class
    predictions = np.argmax(probabilities, axis=1)
    for i in range(len(null_set)):
        s = predictions[i]
        e_id = null_set[i]
        if e_id not in labels:
            labels[e_id] = [sample_keys[s]]
            """
            for assoc in assoc_dict[sample_keys[s]]:
                labels[e_id] += [assoc]
            """
        else:
            labels[e_id].append(sample_keys[s])

    return labels
                    
