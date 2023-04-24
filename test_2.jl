using Plots
using Optim

# 2D Gaussian function
function gaussian_2d(x, y, μx, μy, σx, σy)
    return (1 / (2 * π * σx * σy)) * exp(-((x - μx)^2 / (2 * σx^2) + (y - μy)^2 / (2 * σy^2)))
end

# Maximum likelihood estimation function
function mle(params, photons, σx, σy)
    μx1, μy1, μx2, μy2 = params
    likelihood = 0.0
    for (x, y) in eachrow(photons)
        p1 = gaussian_2d(x, y, μx1, μy1, σx, σy)
        p2 = gaussian_2d(x, y, μx2, μy2, σx, σy)
        likelihood += log(p1 + p2)
    end
    return -likelihood
end

# Parameters for point sources
σx, σy = 3.0, 3.0

# Detected photon positions
photons = [5.2 7.1;
           #4.8 6.9;
           10.1 11.9;
           9.9 12.3;
           #5.3 7.2;
           #10.3 11.8
           ]

# Optimize MLE function to find the most likely positions of the point sources
initial_guess = [4.0, 6.0, 9.0, 11.0]
result = optimize(p -> mle(p, photons, σx, σy), initial_guess, BFGS())
estimated_positions = result.minimizer

# Define the grid for the probability map
x = range(0.0, stop=15.0, length=200)
y = range(0.0, stop=20.0, length=200)

# Calculate the probability distribution given the detected photon positions
probability_map = zeros(length(y), length(x))
for i in 1:length(y)
    for j in 1:length(x)
        p1 = gaussian_2d(x[j], y[i], estimated_positions[1], estimated_positions[2], σx, σy)
        p2 = gaussian_2d(x[j], y[i], estimated_positions[3], estimated_positions[4], σx, σy)
        probability_map[i, j] = p1 + p2
    end
end

# Normalize the probability map
probability_map /= sum(probability_map)

# Plot the probability map
heatmap(x, y, probability_map, title="Probability Map", xlabel="X", ylabel="Y", color=:viridis)
scatter!(photons[:, 1], photons[:, 2], marker=:circle, markersize=8, label="Detected Photons")
scatter!([estimated_positions[1], estimated_positions[3]], [estimated_positions[2], estimated_positions[4]], marker=:star, markersize=10, label="Estimated Positions")
