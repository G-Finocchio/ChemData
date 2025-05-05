# Function that is consistent to most papers
r2 <- function(actual, predicted) {
  mean_actual <- mean(actual)
  
  total_sum_squares <- sum((actual - mean_actual) ^ 2)
  
  residual_sum_squares <- sum((actual - predicted) ^ 2)
  
  r_squared <- 1 - (residual_sum_squares / total_sum_squares)
  #r_squared=corr(actual, predicted)^2 #other option of r2 less used in papers
  
  return(r_squared)
}



element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}

# Function that create the plots of Predicted vs Observed Yields 
plot_predicted_yield<-function(PredictedYield, ObservedYield, title, fn_save)
{
# Calculate R² and RMSE
rsq <- r2(ObservedYield, PredictedYield)
rmse <- sqrt(mean((ObservedYield - PredictedYield)^2))

# Create dataframe
df <- data.frame(ObservedYield, PredictedYield)



# Plot

p<-ggplot(data = df, aes(x = PredictedYield, y = ObservedYield)) + 
  theme(text = element_text(size = 40)) +
  labs(title = title) +
  theme(plot.title = element_textbox(hjust = 0.5, margin = ggplot2::margin(t = 5, b = 5))) +
  geom_point(color = "cornflowerblue", size = 3) +
  xlab("Predicted Yield") + ylab("Observed Yield") +
  geom_abline(intercept = 0, slope = 1, size = 1.5, linetype = "dashed") +
  xlim(-0.5, 1) +  
  annotate("text", x = c(-0.35, -0.28), y = c(0.95, 0.89),
           label = c(paste("R^2:", round(rsq, 6)),
                     paste("RMSE:", round(rmse, 6))),
           size = 10, parse = TRUE)
pdf(fn_save, width=10,height=10)
print(p)
dev.off()
print("The plot is saved in folder Results.")
}




