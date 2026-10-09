### Monte Carlo Simulation and Importance Sampling
import numpy as np

# Parameters for normal distributions
mu1, sigma1 = 10, 2  # Distribution of X1
mu2, sigma2 = 5, 1    # Distribution of X2

# Number of samples
num_samples = 100000

# Generate random samples from normal distributions
X1_samples = np.random.normal(mu1, sigma1, num_samples)
X2_samples = np.random.normal(mu2, sigma2, num_samples)

# Limit state function
def limit_state_function(x1, x2):
    return x1 - x2

# Evaluate the limit state function for all samples
g_values = limit_state_function(X1_samples, X2_samples)

# Calculate failure probability (g < 0)
failure_probability = np.mean(g_values < 0)

print(f"Failure Probability (Simple Monte Carlo): {failure_probability:.6f}")

# Improved Monte Carlo (Importance Sampling)
# New distributions for X1 and X2 (focused on the failure region)
mu1_IS, sigma1_IS = 5, 2  # New distribution for X1
mu2_IS, sigma2_IS = 10, 1  # New distribution for X2

# Generate random samples from the new distributions
X1_samples_IS = np.random.normal(mu1_IS, sigma1_IS, num_samples)
X2_samples_IS = np.random.normal(mu2_IS, sigma2_IS, num_samples)

# Evaluate the limit state function for the new samples
g_values_IS = limit_state_function(X1_samples_IS, X2_samples_IS)

# Calculate Importance Sampling weights
weights = (np.exp(-0.5 * ((X1_samples_IS - mu1) / sigma1)**2) / np.exp(-0.5 * ((X1_samples_IS - mu1_IS) / sigma1_IS)**2) * (np.exp(-0.5 * ((X2_samples_IS - mu2) / sigma2)**2) / np.exp(-0.5 * ((X2_samples_IS - mu2_IS) / sigma2_IS)**2)))

# Calculate failure probability with Importance Sampling
failure_probability_IS = np.mean((g_values_IS < 0) * weights)

print(f"Failure Probability (Improved Monte Carlo): {failure_probability_IS:.6f}")

#-----------------------------------------------------------------------
### Neural Network for Failure Probability Estimation

import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

# 1. Generate training data
num_samples = 10000

# Generate random variables
X1 = np.random.normal(10, 2, num_samples)
X2 = np.random.normal(5, 1, num_samples)

# Calculate the limit state function and labels (0 = failure, 1 = safe)
g = X1 - X2
labels = (g >= 0).astype(int)  # 1 if safe, 0 if failure

# Combine data
data = np.column_stack((X1, X2))

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.2, random_state=42)

# 2. Build the neural network
model = Sequential([
    Dense(16, activation='relu', input_shape=(2,)),  # Input layer with 2 neurons (X1, X2)
    Dense(8, activation='relu'),
    Dense(1, activation='sigmoid')  # Output between 0 and 1 (failure probability)
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# Train the model
history = model.fit(X_train, y_train, epochs=50, batch_size=32, validation_split=0.2)

# 3. Predict failure probability using the neural network
num_mc_samples = 100000
X1_mc = np.random.normal(10, 2, num_mc_samples)
X2_mc = np.random.normal(5, 1, num_mc_samples)
data_mc = np.column_stack((X1_mc, X2_mc))

# Predict failure probability
predictions = model.predict(data_mc)
failure_probability_ann = np.mean(predictions < 0.5)  # If predicted probability < 0.5, it's a failure

# Calculate true failure probability (for comparison)
g_mc = X1_mc - X2_mc
failure_probability_true = np.mean(g_mc < 0)

print(f"True Failure Probability (Monte Carlo): {failure_probability_true:.6f}")
print(f"Predicted Failure Probability (Neural Network): {failure_probability_ann:.6f}")

# Plot training accuracy
plt.plot(history.history['accuracy'], label='Train Accuracy')
plt.plot(history.history['val_accuracy'], label='Validation Accuracy')
plt.xlabel('Epoch')
plt.ylabel('Accuracy')
plt.legend()
plt.show()

#-----------------------------------------------------------------------
### Response Surface Method (RSM) with Polynomial Regression

import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

# 1. Generate training data
num_samples = 1000

# Generate random variables
X1 = np.random.normal(2, 0.5, num_samples)
X2 = np.random.normal(3, 0.8, num_samples)

# True limit state function
def true_limit_state(x1, x2):
    return x1**3 - 2*x2**2 + 5

# Evaluate the true limit state function
g_true = true_limit_state(X1, X2)

# 2. Build the Response Surface Model (quadratic polynomial)
# Create polynomial features
poly = PolynomialFeatures(degree=2)
X_data = np.column_stack((X1, X2))
X_poly = poly.fit_transform(X_data)

# Train the regression model
model = LinearRegression()
model.fit(X_poly, g_true)

# Response surface function
def response_surface(x1, x2):
    X_input = poly.transform(np.column_stack((x1, x2)))
    return model.predict(X_input)

# 3. Calculate failure probability using the response surface model
num_mc = 100000
X1_mc = np.random.normal(2, 0.5, num_mc)
X2_mc = np.random.normal(3, 0.8, num_mc)

# Predict using the response surface model
g_pred = response_surface(X1_mc, X2_mc)
failure_prob_rsm = np.mean(g_pred < 0)

# Calculate true failure probability (Monte Carlo)
g_true_mc = true_limit_state(X1_mc, X2_mc)
failure_prob_true = np.mean(g_true_mc < 0)

print(f"True Failure Probability (Monte Carlo): {failure_prob_true:.6f}")
print(f"Predicted Failure Probability (Response Surface): {failure_prob_rsm:.6f}")

# 4. Visualize the response surface
x1_range = np.linspace(0, 4, 100)
x2_range = np.linspace(1, 5, 100)
X1_grid, X2_grid = np.meshgrid(x1_range, x2_range)
g_grid = response_surface(X1_grid.ravel(), X2_grid.ravel()).reshape(X1_grid.shape)

plt.figure(1, figsize=(20, 10))
plt.contourf(X1_grid, X2_grid, g_grid, levels=20, cmap='viridis')
plt.colorbar(label='Limit State Function (g)')
plt.scatter(X1_mc[g_true_mc < 0], X2_mc[g_true_mc < 0], c='red', s=1, alpha=0.3, label='Failure')
plt.scatter(X1_mc[g_true_mc >= 0], X2_mc[g_true_mc >= 0], c='blue', s=1, alpha=0.3, label='Safe')
plt.xlabel('X1')
plt.ylabel('X2')
plt.title('Response Surface and Data Distribution')
plt.legend()
plt.show()