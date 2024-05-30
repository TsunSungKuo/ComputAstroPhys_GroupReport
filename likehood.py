import numpy as np
import matplotlib.pyplot as plt

def rbf_kernel(X1, X2, length_scale, variance):
    sqdist = np.sum(X1**2, 1).reshape(-1, 1) + np.sum(X2**2, 1) - 2 * np.dot(X1, X2.T)
    return variance * np.exp(-0.5 / length_scale**2 * sqdist)

def log_marginal_likelihood(theta, X, y):
    length_scale, variance, noise_variance = theta
    K = rbf_kernel(X, X, length_scale, variance) + noise_variance * np.eye(len(X))
    L = np.linalg.cholesky(K)
    alpha = np.linalg.solve(L.T, np.linalg.solve(L, y))
    log_likelihood = -0.5 * np.dot(y.T, alpha)
    log_likelihood -= np.sum(np.log(np.diag(L)))
    log_likelihood -= len(X) / 2 * np.log(2 * np.pi)
    return log_likelihood

# Example data for a shock tube problem
# Spatial positions along the shock tube
X = np.array([[0.1], [0.2], [0.3], [0.4], [0.5], [0.6], [0.7], [0.8], [0.9]])

# Corresponding pressure values at those positions
y = np.array([1.0, 0.9, 0.7, 0.4, 0.2, 0.1, 0.05, 0.02, 0.01])

# Initial hyperparameters
theta = [0.1, 1.0, 1e-8]

# Compute the log marginal likelihood
log_likelihood = log_marginal_likelihood(theta, X, y)
print("Log Marginal Likelihood:", log_likelihood)

# Fit and predict using manual GP
K = rbf_kernel(X, X, length_scale=theta[0], variance=theta[1]) + theta[2] * np.eye(len(X))
L = np.linalg.cholesky(K)
alpha = np.linalg.solve(L.T, np.linalg.solve(L, y))

# Test data: spatial positions where we want to predict pressure
X_test = np.linspace(0, 1, 100).reshape(-1, 1)
K_s = rbf_kernel(X, X_test, length_scale=theta[0], variance=theta[1])
K_ss = rbf_kernel(X_test, X_test, length_scale=theta[0], variance=theta[1]) + theta[2] * np.eye(len(X_test))

# Mean prediction
mu_s = K_s.T.dot(alpha)

# Covariance of the prediction
v = np.linalg.solve(L, K_s)
cov_s = K_ss - v.T.dot(v)
sigma_s = np.sqrt(np.diag(cov_s))

# Plot the results
plt.figure(figsize=(10, 5))
plt.plot(X, y, 'r.', markersize=10, label='Observed data')
plt.plot(X_test, mu_s, 'b-', label='GP prediction')
plt.fill_between(X_test.ravel(), mu_s - 1.96 * sigma_s, mu_s + 1.96 * sigma_s, alpha=0.2, color='blue', label='95% confidence interval')
plt.xlabel('Position along the shock tube')
plt.ylabel('Pressure')
plt.legend()
plt.show()
