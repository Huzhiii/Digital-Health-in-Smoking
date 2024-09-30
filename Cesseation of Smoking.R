# Load required libraries
library(gemtc)       # For Bayesian network meta-analysis
library(netmeta)     # For frequentist network meta-analysis
library(meta)        # For conventional meta-analysis
library(metafor)     # For advanced meta-analysis models
library(ggplot2)     # For data visualization
library(readxl)      # For reading Excel files

# Import data
data <- read_excel("频率point.xlsx")


# Convert columns to numeric data types
data$responders <- as.numeric(data$responders)
data$sampleSize <- as.numeric(data$sampleSize)

# Pairwise comparisons for the treatments
p1 <- pairwise(treatment, event = responders, n = sampleSize,
               studlab = study, data = data, sm = "RR")

# Frequentist Network Meta-analysis using the random effects model
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data = p1, sm = "RR",
                comb.fixed = FALSE, comb.random = TRUE, reference = "Standard_Care")

# Summary of network meta-analysis
summary(net1)

# Create a forest plot
forest(net1, ref = "Standard_Care",
       pooled = "random", digits = 2,
       smlab = "Random effects model",
       xlab = "Relative Risk",
       leftlabs = "Digital Intervention")

# Plot Bayesian network diagram
network <- mtc.network(data)
plot(network, use.description = TRUE, vertex.label.cex = 1,
     vertex.size = 5, vertex.shape = "circle", vertex.label.color = "darkblue",
     vertex.label.dist = 1, vertex.label.degree = -pi/2, vertex.label.cex = 200,
     vertex.color = "red",
     dynamic.edge.width = TRUE,
     edge.color = "blue",
     vertex.label.font = 2)


# Create a league table and export as CSV
netleague <- netleague(net1, bracket = "(", digits = 2)
write.csv(netleague$random, "netleague.csv")

# Generate P-scores for treatment ranking
netrank(net1, small.values = "bad")

# Generate heatmap to assess inconsistency
netheat(net1, random = TRUE)

# Display network consistency check
netsplit(net1)
sink("output.txt")          # Redirect output to a file
print(netsplit(net1))       # Save network split details
sink()

# Publication Bias Analysis using Funnel Plot
colors <- rainbow(45)      # Generate 45 unique colors
pch_values <- rep(0:25, length.out = 45)  # Set point character values

# Create funnel plot with treatment order
funnel(net1,
       order = c("Placebo", "Standard_Care", "Face_to_Face", "Phone", "SMS",
                 "Email", "Web", "App", "Multicomponent_Intervention", "Interactive_SMS",
                 "Interactive_Web", "Interactive_app", "Customized_SMS", "Customized_Email",
                 "Customized_Web", "Customized_App", "Group_Customized_SMS",
                 "Group_Customized_App", "Group_Customized_Web", "Group_Customized_Phone",
                 "Group_Customized_Multicomponent_Intervention"),
       pch = pch_values,
       col = colors[1:45],  # Use 45 unique colors
       linreg = TRUE,
       legend = FALSE)

# Further Funnel Plots for different treatment comparisons (as per study requirements)


# Regression Analysis for various covariates using Bayesian network meta-analysis
data <- read_excel("Frequency_Calculations.xlsx")

data$responders <- as.numeric(data$responders)
data$sampleSize <- as.numeric(data$sampleSize)

network <- mtc.network(data)
plot(network, use.description = TRUE, vertex.label.cex = 1,
     vertex.size = 5, vertex.shape = "circle", vertex.label.color = "darkblue",
     vertex.label.dist = 1, vertex.label.degree = -pi/2, vertex.label.cex = 200,
     vertex.color = "red",
     dynamic.edge.width = TRUE,
     edge.color = "blue",
     vertex.label.font = 2)

# Bayesian modeling and regression analysis
model <- mtc.model(network, type = "consistency", n.chain = 4, likelihood = "binom", link = "log", linearModel = "random")

results <- mtc.run(model, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results)
forest(relative.effect(results, "Standard_Care"))

# Rank probability plot
rank <- rank.probability(results, preferredDirection = -1)
plot(rank, beside = TRUE)
forest(results, t1 = "Standard_Care", use.description = TRUE)

# Analysis for other variables (sex, age, smoking status, year, etc.)
studies <- data

studies$sm <- as.numeric(studies$sm)
studies$Sex <- as.numeric(studies$Sex)
studies <- studies[!is.na(studies$Average_age), ]
studies$`2010` <- ifelse(tolower(studies$`2010`) == "yes", 1, ifelse(tolower(studies$`2010`) == "no", 0, NA))
studies$`2015` <- ifelse(tolower(studies$`2015`) == "yes", 1, ifelse(tolower(studies$`2015`) == "no", 0, NA))
studies$bioreport <- ifelse(tolower(studies$bioreport) == "yes", 1, ifelse(tolower(studies$bioreport) == "no", 0, NA))
studies$drug <- ifelse(tolower(studies$drug) == "yes", 1, ifelse(tolower(studies$drug) == "no", 0, NA))
studies$money <- ifelse(tolower(studies$money) == "yes", 1, ifelse(tolower(studies$money) == "no", 0, NA))

network <- mtc.network(data = studies, studies = studies)


# Running regression models for different variables

# Model for 'Drug' variable
model_drug <- mtc.model(network, type = "regression",
                        regressor = list(coefficient = 'shared', 
                                         variable = 'drug', 
                                         control = 'Standard_Care'))
results_drug <- mtc.run(model_drug, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results_drug)

# Model for 'Time' variable (assuming 'Time' as another variable in your dataset)
model_time <- mtc.model(network, type = "regression", 
                        regressor = list(coefficient = 'shared', 
                                         variable = 'Time', 
                                         control = 'Standard_Care'))
results_time <- mtc.run(model_time, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results_time)

# Model for 'Money' variable
model_money <- mtc.model(network, type = "regression", 
                         regressor = list(coefficient = 'shared', 
                                          variable = 'money', 
                                          control = 'Standard_Care'))
results_money <- mtc.run(model_money, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results_money)

summary(results_country)


# Model for 'Bio-reporting' variable
model_bio <- mtc.model(network, type = "regression", 
                       regressor = list(coefficient = 'shared', 
                                        variable = 'bioreport', 
                                        control = 'Standard_Care'))
results_bio <- mtc.run(model_bio, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results_bio)

# Model for 'Sex' variable
model_sex <- mtc.model(network, type = "regression", 
                       regressor = list(coefficient = 'shared', 
                                        variable = 'Sex', 
                                        control = 'Standard_Care'))
results_sex <- mtc.run(model_sex, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results_sex)

# Model for 'Average Age'
model_age <- mtc.model(network, type = "regression", 
                       regressor = list(coefficient = 'shared', 
                                        variable = 'Average_age', 
                                        control = 'Standard_Care'))
results_age <- mtc.run(model_age, n.adapt = 5000, n.iter = 20000, thin = 1)
summary(results_age)

# Model for 'Smoking Status'
model_sm <- mtc.model(network, type = "regression", 
                      regressor = list(coefficient = 'shared', 
                                       variable = 'sm', 
                                       control = 'Standard_Care'))
results_sm <- mtc.run(model_sm, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results_sm)

# Model for studies conducted after 2010
model_2010 <- mtc.model(network, type = "regression", 
                        regressor = list(coefficient = 'shared', 
                                         variable = '2010', 
                                         control = 'Standard_Care'))
results_2010 <- mtc.run(model_2010, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results_2010)

# Model for studies conducted after 2015
model_2015 <- mtc.model(network, type = "regression", 
                        regressor = list(coefficient = 'shared', 
                                         variable = '2015', 
                                         control = 'Standard_Care'))
results_2015 <- mtc.run(model_2015, n.adapt = 20000, n.iter = 50000, thin = 1)
summary(results_2015)

# Final Step: Summarizing the results from all models
# The following lines save all summaries to an output file for documentation purposes

sink("regression_results.txt")          # Redirect output to a file
print("Summary for Drug Model:")
print(summary(results_drug))

print("Summary for Time Model:")
print(summary(results_time))

print("Summary for Money Model:")
print(summary(results_money))

print("Summary for Bio-reporting Model:")
print(summary(results_bio))

print("Summary for Sex Model:")
print(summary(results_sex))

print("Summary for Age Model:")
print(summary(results_age))

print("Summary for Smoking Model:")
print(summary(results_sm))

print("Summary for 2010 Model:")
print(summary(results_2010))

print("Summary for 2015 Model:")
print(summary(results_2015))
sink()    # End of redirection, reverting to normal console output


