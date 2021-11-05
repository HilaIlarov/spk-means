import spkmeans
import numpy as np
import pandas as pd
import sys
from numpy import linalg as LA

def spk(k, goal, file_name):
    max_iter = 300
    np.random.seed(0)
    #Open file
    try:
        file = open(file_name, "r")
    except:
        print("An Error Has Occured")
        return None
    #takes first line of the file and checks how many words in it
    d = len(file.readline().split(','))   
    lst_mu = []
    #Create data frame
    data = pd.read_csv(file_name, header = None)
    N = data.index[-1] + 1
    lst_vectors = data.values.tolist()
    # calling C by fit
    datapoints = spkmeans.send_T_to_python(k, goal, N, d, lst_vectors)
    #now we have initial coordinates. We need to send it to C and activate kMeans.
    #running kmeans++
    N = len(datapoints)
    k = len(datapoints[0])
    observations = []   #observations = list of random indexes for initial centroids
    rand_index = np.random.choice(N,1)[0]
    observations.append(int(rand_index))
    vector_mu = datapoints[rand_index]
    lst_mu.append(vector_mu)
    D = [0 for i in range(N)]   #D = list of min values
    P =[0 for i in range(N)]    #P = list of probabilities
    Z = 1
    while (Z != k):
        #for each vector x_i (for 1<=i<=N) we compute the norm with mu_j for 1<=j<=Z
        for i in range(N):
            dp_i = datapoints[i]
            min_val = min(LA.norm(np.array(dp_i) - np.array(mu)) ** 2 for mu in lst_mu)
            D[i] = min_val
        Z += 1
        #calculating probabilities
        sum_D = sum(D)
        for i in range(N):
            P[i] = D[i]/(sum_D)
        #choose the next random key by the probabilities
        rand_index = np.random.choice(N, 1, p = P)[0]
        observations.append(int(rand_index))
        vector_mu = datapoints[rand_index]
        lst_mu.append(vector_mu)
    d = k     #calling C by fit (HW2)
    centroids = spkmeans.fit(N, d, k, max_iter, datapoints, lst_mu)
    print(','.join(str(x) for x in observations))
    centroids = np.array(centroids)
    #rounding the results
    for i in range(len(centroids)):
        centroid = np.round(centroids[i], 4)
        centroids[i] = centroid
    centroids = centroids.tolist()
    #printing the centroids
    for cen in centroids:
        for i in range(len(cen)):
            value = cen[i] + 0.0
            if(i < len(cen) - 1):
                print("{:.4f}".format(value), end = ",")
            else:
                print("{:.4f}".format(value), end = "")
        print()

def main():
    k_str = sys.argv[1]
    k = int(k_str)
    goal = sys.argv[2]  
    file_name = sys.argv[3]
    spk(k, goal, file_name)

main()