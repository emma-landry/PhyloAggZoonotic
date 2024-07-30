library(MASS)
comm_subset <- comm[, intersect(colnames(comm), vZoo$virus), drop = FALSE]
comm_subset <- comm_subset[rowSums(comm_subset) != 0, ]

phydist.mini.subset <- phydist.mini[intersect(rownames(comm_subset), colnames(phydist.mini)),
                                    intersect(rownames(comm_subset), colnames(phydist.mini)),
                                    drop = FALSE]

k_vals <- 2:10
stress_vals <- c()
for (i in k_vals) {
  fit <- isoMDS(phydist.mini.subset, k = i)
  stress_vals <- c(stress_vals, fit$stress)
}

plot(k_vals, stress_vals, type = "l")

better <- c()
for (i in k_vals) {
  fit <- isoMDS(phydist.mini.subset, k = i)
  
  reduced <- as.matrix(fit$points)
  
  predictors <- t(comm_subset) %*% reduced
  
  vZoo_subset <- vZoo[vZoo$virus %in% intersect(colnames(comm), vZoo$virus), ]
  outcome <- vZoo_subset$zoo
  
  outcome <- as.factor(outcome)
  data_for_glm <- data.frame(outcome, predictors)
  logistic_model <- glm(outcome ~ ., data = data_for_glm, family = binomial)
  
  summary(logistic_model)  
  
  comparison1 <- data.frame(MDS = logistic_model$fitted.values)
  comparison1 <- comparison1[rownames(comparison1) %in% intersect(rownames(comparison1), z_$virus), , drop = F]
  comparison2 <- data.frame(MDS = comparison1$MDS,
                            other = model.9$fitted.values)
  rownames(comparison2) <- rownames(comparison1)  
  comparison2$true <- z_$zoo
  comparison2$MDSBetter <- abs(comparison2$MDS - comparison2$true) < abs(comparison2$other - comparison2$true)
  better <- c(better, sum(comparison2$MDSBetter))
  
}

hosts_virus8 <- which(comm_subset[, 8] == 1)
distances8 <- phydist.mini.subset[hosts_virus8, hosts_virus8]

mds8 <- isoMDS(distances8, k = 2)

plot(mds8$points)

### bar plots ----------------------------------
unique_orders <- unique(v$vOrder)
unique_orders <- unique_orders[order(unique_orders)]
unique_family <- unique(v$vFamily)
unique_family <- unique_family[order(unique_family)]
unique_subfamily <- unique(v$vSubfamily)
unique_subfamily <- unique_subfamily[order(unique_subfamily)]

orders_distances <- c(0.0206, 0.0205, 0.0093, 0.0285, 0.0205)
orders_vals <- data.frame("Order" = unique_orders,
                          "Distance" = orders_distances,
                          "Zoonotic" = rep(0, length(orders_distances)))

i <- 0
for (j in unique_orders) {
  i <- i + 1
  subdata <- v[v$vOrder == j, ]
  prop <- sum(subdata$IsZoonotic) / nrow(subdata)
  orders_vals$Zoonotic[i] <- prop
}

family_distances <- c(0.0146, 0.0333, 0.0251, 0.0177, 0.0156, 0.0015,
                      0.0095, 0.0203, 0.0071, 0.0000, 0.0074, 0.0150,
                      0.0202, 0.0327, 0.0170, 0.0215, 0.0142, 0.0093,
                      0.0260, 0.0224, 0.0117, 0.0285, 0.0290, 0.0199,
                      0.0195, 0.0228, 0.0174, 0.0151, 0.0084)
family_vals <- data.frame("Family" = unique_family,
                          "Distance" = family_distances,
                          "Zoonotic" = rep(0, length(family_distances)))

i <- 0
for (j in unique_family) {
  i <- i + 1
  subdata <- v[v$vFamily == j, ]
  prop <- sum(subdata$IsZoonotic) / nrow(subdata)
  family_vals$Zoonotic[i] <- prop
}

subfamily_distances <- c(0.0183, 0.0182, 0.0199, 0.0069, 0.0198,
                         0.0226, 0.0277, 0.0145, 0.0142, 0.0194,
                         0.0236, 0.0224, 0.000)
subfamily_vals <- data.frame("Subfamily" = unique_subfamily[-14],
                          "Distance" = subfamily_distances,
                          "Zoonotic" = rep(0, length(subfamily_distances)))

i <- 0
for (j in unique_subfamily[-14]) {
  i <- i + 1
  notNA <- v[!is.na(v$vSubfamily), ]
  subdata <- notNA[notNA$vSubfamily == j, ]
  prop <- sum(subdata$IsZoonotic) / nrow(subdata)
  subfamily_vals$Zoonotic[i] <- prop
}

ggplot(orders_vals, aes(x = Distance, y = reorder(Order, Distance), fill = Zoonotic)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    title = "Cosine distance across virus orders",
    x = "Cosine Distance",
    y = "Order",
    fill = "Prop. Zoonotic"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    legend.position = c(0.95, 0.45),  # Position the legend inside the plot at the bottom right
    legend.justification = c("right", "bottom")  # Adjust legend justification
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("cosine_distance_across_virus_orders.pdf", width = 10, height = 7)

ggplot(family_vals, aes(x = Distance, y = reorder(Family, Distance), fill = Zoonotic)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    title = "Cosine distance across virus families",
    x = "Cosine Distance",
    y = "Family",
    fill = "Prop. Zoonotic"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    legend.position = c(0.95, 0.45),  # Position the legend inside the plot at the bottom right
    legend.justification = c("right", "bottom")  # Adjust legend justification
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("cosine_distance_across_family_orders.pdf", width = 10, height = 7)

ggplot(subfamily_vals, aes(x = Distance, y = reorder(Subfamily, Distance), fill = Zoonotic)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    title = "Cosine distance across virus subfamilies",
    x = "Cosine Distance",
    y = "Subfamily",
    fill = "Prop. Zoonotic"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),  # Center the plot title
    legend.position = c(0.95, 0.45),  # Position the legend inside the plot at the bottom right
    legend.justification = c("right", "bottom")  # Adjust legend justification
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("cosine_distance_across_subfamily_orders.pdf", width = 10, height = 7)