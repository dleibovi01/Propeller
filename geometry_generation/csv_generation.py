import numpy as np



radiusq = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1])


# beta=np.array([24.65, 25.65, 26.35, 26.75, 26.49, 25.57, 23.95, 21.72, 19.38, 17.53, 15.8, 14.4, 13.19, 12.09, 11.38, 10.13, 8.7, 6.5, 4.26])
Cq = np.array([0.157, 0.169, 0.181, 0.183, 0.183, 0.182, 0.18, 0.176, 0.173, 0.17, 0.167, 0.162, 0.156, 0.15, 0.142, 0.132, 0.117, 0.095, 0.073])


beta0 = np.array([33.52, 37.73, 39.74, 36.96, 33.71, 30.82, 28.05, 25.25, 22.79, 20.71, 18.88, 17.24, 15.79, 14.54, 13.40, 12.30, 11.40, 10.33, 9.26]) # APC Slow Flyer 9x6
beta1 = np.array([31.0, 35.61, 39.44, 42.45, 41.18, 37.35, 33.74, 30.77, 27.98, 25.54, 23.41, 21.87, 20.17, 18.87, 17.76, 16.05, 14.17, 12.29, 10.41]) # APC Sport propeller 9x7
beta2 = np.array([30.15, 27.15, 20.45, 20.95, 21.26, 21.60, 21.94, 22.14, 22.14, 21.78, 21.22, 20.29, 19.18, 17.87, 16.65, 15.64, 14.52, 12.8, 11.04]) # GWS Slow Flyer
beta3 = np.array([31.5, 30.52, 29.64, 27.42, 25.57, 24.0, 22.55, 21.25, 20.03, 18.75, 17.57, 16.35, 15.1, 13.81, 12.49, 11.39, 10.14, 8.5, 6.84])

beta_avg = np.average(np.concatenate((beta0.reshape(-1, 1), beta1.reshape(-1, 1), beta2.reshape(-1, 1), beta3.reshape(-1, 1)), axis=1), axis=1)

# print(beta_avg.shape)
# print(beta_avg)
# input()
# print(beta0.shape)
# print(beta1.shape)
# print(beta2.shape)
# print(beta3.shape)

skew0 = np.zeros(19)
skew1 = np.array([-1.8, -2.7, -3.32, -3.99, -4.39, -4.39, -4.39, -3.65, -3.14, -1.98, -0.82, 0.81, 2.49, 4.47, 6.33, 8.48, 10.64, 12.93, 14.1]) # KVLCC KP458
skew2 = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3.0, 2.5, 2, 1.5, 1, 0.5, 0, -0.5, -1, -1.5, -2.0, -2.5, -3.0])

skew_avg = np.average(np.concatenate((skew0.reshape(-1, 1), skew1.reshape(-1, 1), skew2.reshape(-1, 1)), axis=1), axis=1)

c_beta = np.random.rand(4, 1)
noise = np.random.normal(0, 1, size = beta0.shape)
# c_beta = np.array([0, 0, 1, 0])
beta = (c_beta[0] * beta0 + c_beta[1] * beta1 + c_beta[2] * beta2 + c_beta[3] * beta3) / np.sum(c_beta)  + 0.05 * noise * beta_avg


c_skew = np.random.rand(3, 1)
noise = np.random.normal(0, 1, size = beta0.shape)
# c_skew = np.array([1, 0, 0])
Skewsq = (c_skew[0] * skew0 + c_skew[1] * skew1 + c_skew[2] * skew2) / np.sum(c_skew) + 0.05 * noise * skew_avg

# print(c_beta)
# print(beta)


# print(Skewsq)

# Skewsq = np.zeros(19)
NACA1s = np.zeros(19)
NACA2s = np.zeros(19)


data = np.stack((radiusq, beta, Cq, Skewsq, NACA1s, NACA2s), axis=1)

# data = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
header = "radiusq,beta,Cq,Skewsq,NACA1s,NACA2s"

np.savetxt("case_info.csv", data, delimiter=",", header=header, comments="")