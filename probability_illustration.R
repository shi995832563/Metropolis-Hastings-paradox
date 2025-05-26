par(mfrow=c(1,1),mar=c(3,3.6,1.5,0),oma=c(0,0,0,1),mgp=c(2,0.8,0)) 

x_values <- seq(-5, 5, length.out = 100)
v2 <- 0.8
# Calculate the density for each distribution
density_standard <- dnorm(x_values, mean = 0, sd = 1)
density_normal_1 <- dnorm(x_values, mean = 1, sd = sqrt(v2))

# Create the plot for the standard normal distribution
plot(x_values, density_standard, type = "l", col = "red", lwd = 2,
     xlab = "", ylab = "Probability density function", ylim=c(0,0.5),
     main = "")


# Add a vertical line at the mean of the second distribution
#abline(v = 1, col = "blue", lty = 2)
#abline(v = -1, col = "blue", lty = 2)
#abline(v = 0, col = "red", lty = 2)
#segments(x0 = 1, y0 = 0, y1 = dnorm(1, mean = 1, sd = sqrt(v2)), col = "blue", lty = 2, lwd = 2)
#segments(x0 = -1, y0 = 0, y1 = dnorm(-1, mean = 1, sd = sqrt(v2)), col = "blue", lty = 2, lwd = 2)


x_fill <- seq(-1, 1, length.out = 100)
y_fill <- dnorm(x_fill, mean = 1, sd = sqrt(v2))
polygon(c(-1, x_fill, 1), c(0, y_fill, 0), col = "lightblue", border = NA)

x_fill <- seq(-5, -1, length.out = 100)
y_fill <- dnorm(x_fill, mean = 1, sd = sqrt(v2))
polygon(c(-5, x_fill, -1), c(0, y_fill, 0), col = "pink", border = NA)

x_fill <- seq(1, 5, length.out = 100)
y_fill <- dnorm(x_fill, mean = 1, sd = sqrt(v2))
polygon(c(1, x_fill, 5), c(0, y_fill, 0), col = "pink", border = NA)

# Annotate the mean
text(1.3, -0.03, expression(math(theta)[i]^{(t-1)}), col = "blue", pos = 3)
#text(-0.8, -0.03, expression(-x[i]^{(t-1)}), col = "blue", pos = 3)

# Add notation "I" in the middle of the light blue area
text(0.3, 0.17, "H", col = "black", cex = 1.5)  # Adjust y position as needed
text(1.5, 0.17, "L", col = "black", cex = 1.5)  # Adjust y position as needed

# Add notation "I" in the middle of the light blue area
text(3, 0.47, "Proposal distribution", col = "blue", cex = 1.5)  # Adjust y position as needed
text(-2.5, 0.4, "Evaluation distribution", col = "red", cex = 1.5)  # Adjust y position as needed

lines(x_values, density_normal_1, col = "blue", lwd = 2)
lines(x_values, density_standard, type = "l", col = "red", lwd = 2)

segments(x0 = 1, y0 = 0, y1 = dnorm(1, mean = 0, sd = sqrt(v1)), col = "black", lty = 2, lwd = 2)
segments(x0 = -1, y0 = 0, y1 = dnorm(1, mean = 0, sd = sqrt(v1)), col = "black", lty = 2, lwd = 2)
segments(x0=-1, y0=max(dnorm(1, mean = 0, sd = sqrt(v1))),x1=1, y1=max(dnorm(1, mean = 0, sd = sqrt(v1))), col = "black", lty = 2, lwd = 2)
