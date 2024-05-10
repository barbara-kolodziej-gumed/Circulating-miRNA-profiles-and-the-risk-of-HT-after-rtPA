

################################################################################

#### DATASET BEFORE MATCHING

library(MatchIt)
library(readxl)
library(dplyr)
library(writexl)
library(ggplot2)
library(gridExtra)
library(rlang)
options(scipten = 999)

Samples_4_matching <- read_excel("C:/Users/adawys/.../Samples_4_matching.xlsx", 
                                 col_types = c("text", "text", "numeric", 
                                               "text", "numeric", "text", "text"))
Samples_4_matching$Hemorrhage <- factor(Samples_4_matching$Hemorrhage)

#No matching; constructing a pre-match matchit object
m.out0 <- matchit(Hemorrhage ~ Age + Gender + NIHSS + OCSP + Proteomics, 
                  data = Samples_4_matching,
                  method = NULL, distance = "glm")


# 1:1 Nearesr Neigbour PS matching without replacement
NN<- matchit(Hemorrhage ~ Age + Gender + NIHSS + OCSP + Proteomika, 
             data = Samples_4_matching,
             method = "nearest", distance = "glm",
             estimand = "ATT")

# Distribution of Propensity Scores
plot(NN, type = "jitter", interactive = FALSE)

# Balance on the covariates 
plot(NN, type = "density", interactive = FALSE,
     which.xs = ~Age + Gender + NIHSS + OCSP + Proteomika) #treated (black), control (gray)

# Matched data
matched_data <- match.data(NN)


write_xlsx(matched_data, "C:/Users/adawys/.../matched_data.xlsx")





################################################################################


# DATASET AFTER MATCHING WITH MICRORNA VARIABLES
miRNA_non_ich <- read_excel("C:/Users/adawys/Desktop/PROJEKTY/Neurologia/Marcin_Stanczak/miRNA/Baza_danych/miRNA-data.xlsx", 
                         sheet = "non-ich",
                         na="NA")
miRNA_non_ich$Group_level <- rep('Non Ich', 10)

miRNA_ich <- read_excel("C:/Users/adawys/Desktop/PROJEKTY//Neurologia/Marcin_Stanczak/miRNA/Baza_danych/miRNA-data.xlsx", 
                            sheet = "ich",
                        na="NA")
miRNA_ich$Group_level <- rep('Ich', 10)


data <- rbind(miRNA_ich,miRNA_non_ich)

data <- data%>%
  filter(`Gene Name` != "Ratio_47") #exclusion of one observation





# List of variable names
variable_names <- c("`hsa-let-7a-5p`",   "`hsa-miR-1-3p`",    "`hsa-miR-100-5p`",  "`hsa-miR-106b-5p`", "`hsa-miR-10b-5p`",  "`hsa-miR-122-5p`" ,
                    "`hsa-miR-124-3p`",  "`hsa-miR-125b-5p`", "`hsa-miR-126-3p`",  "`hsa-miR-133a-3p`", "`hsa-miR-133b`",    "`hsa-miR-134-5p`",  "`hsa-miR-141-3p`", 
                    "`hsa-miR-143-3p`",  "`hsa-miR-146a-5p`", "`hsa-miR-150-5p`",  "`hsa-miR-155-5p`",  "`hsa-miR-17-5p`",   "`hsa-miR-17-3p`",   "`hsa-miR-18a-5p`", 
                    "`hsa-miR-192-5p`",  "`hsa-miR-195-5p`",  "`hsa-miR-196a-5p`", "`hsa-miR-19a-3p`",  "`hsa-miR-19b-3p`",  "`hsa-miR-200a-3p`", "`hsa-miR-200b-3p`",
                    "`hsa-miR-200c-3p`", "`hsa-miR-203a-3p`", "`hsa-miR-205-5p`",  "`hsa-miR-208a-3p`", "`hsa-miR-20a-5p`",  "`hsa-miR-21-5p`",   "`hsa-miR-210-3p`", 
                    "`hsa-miR-214-3p`",  "`hsa-miR-215-5p`",  "`hsa-miR-221-3p`",  "`hsa-miR-222-3p`",  "`hsa-miR-223-3p`",  "`hsa-miR-224-5p`",  "`hsa-miR-23a-3p`", 
                    "`hsa-miR-25-3p`",   "`hsa-miR-27a-3p`",  "`hsa-miR-296-5p`",  "`hsa-miR-29a-3p`",  "`hsa-miR-30d-5p`",  "`hsa-miR-34a-5p`",  "`hsa-miR-375-3p`", 
                    "`hsa-miR-423-5p`",  "`hsa-miR-499a-5p`", "`hsa-miR-574-3p`",  "`hsa-miR-885-5p`",  "`hsa-miR-9-5p`",    "`hsa-miR-92a-3p`",  "`hsa-miR-93-5p`",  
                    "`hsa-let-7c-5p`",   "`hsa-miR-107`",     "`hsa-miR-10a-5p`",  "`hsa-miR-128-3p`",  "`hsa-miR-130b-3p`", "`hsa-miR-145-5p`",  "`hsa-miR-148a-3p`",
                    "`hsa-miR-15a-5p`",  "`hsa-miR-184`",     "`hsa-miR-193a-5p`", "`hsa-miR-204-5p`",  "`hsa-miR-206`",     "`hsa-miR-211-5p`",  "`hsa-miR-26b-5p`", 
                    "`hsa-miR-30e-5p`",  "`hsa-miR-372-3p`",  "`hsa-miR-373-3p`",  "`hsa-miR-374a-5p`", "`hsa-miR-376c-3p`", "`hsa-miR-7-5p`",    "`hsa-miR-96-5p`",  
                    "`hsa-miR-103a-3p`", "`hsa-miR-15b-5p`",  "`hsa-miR-16-5p`",   "`hsa-miR-191-5p`",  "`hsa-miR-22-3p`",   "`hsa-miR-24-3p`" ,  "`hsa-miR-26a-5p`", 
                    "`hsa-miR-31-5p`")


data$Group_level_new <- ifelse(data$Group_level=="Ich","HT","non-HT")


################################################################################

#### GENERATING PLOTS

clean_label <- function(var_name) {
  # Replace or remove characters as needed
  label <- gsub("`", "", var_name)  # Remove backticks
  return(label)
}


# Loop through each variable
for (i in seq_along(variable_names)) {
  var <- variable_names[i]
  # Clean the variable name for the axis label
  clean_label_var <- clean_label(var)
  
  # Density plot
  plot1 <-ggplot(data, aes_string(x = var, 
                                  fill = "Group_level_new")) +
    geom_density(alpha = 0.7) + # Creates density plot
    labs(title = paste("Density Plot of", clean_label_var), x = clean_label_var, y = "Density", fill="Group") +
    theme_minimal() +
    theme(text = element_text(size = 12), 
          axis.title = element_text(size = 14), 
          plot.title = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))+
    scale_fill_brewer(palette = "Set1") # Adjusts color scheme
  
  
  
  
  # Determine the amount of jitter in the y-axis
  y_jitter_height <- ifelse(all(data[[var]] == 0), 0, 0.1)
  # Create the plot
  plot2 <- ggplot(data, aes_string(x = "Group_level_new", y = var)) +
    geom_boxplot(outlier.shape = NA, width=0.7)+
    geom_jitter(width = 0.1, height = y_jitter_height,  color = "blue", alpha = 0.7, size = 3) +
    labs(title = paste("Boxplot of", clean_label_var), y = clean_label_var, x = "Group")+
    theme_minimal() +
    theme(text = element_text(size = 24), 
          axis.title = element_text(size = 28), 
          plot.title = element_text(size = 28),
          axis.text.x = element_text(size = 24),
          axis.text.y = element_text(size = 24),
          legend.title = element_text(size = 24),
          legend.text = element_text(size = 24))+
    scale_fill_brewer(palette = "Set2")+
    ylim(0,NA)
  
  # Specify the file name
  setwd("C:/Users/adawys/.../")
  file_name <- paste0("Variable_nr_",i, "_",var, "_plot.jpeg")
  
  # Open JPEG device
  jpeg(file_name ,width =6, height = 6, units = "in", res = 1000, pointsize = 35)
  
  # Print the plot
  combined_plot <- grid.arrange(plot2)
  
  # Close the device
  dev.off()
}

################################################################################

#### GENERATING A TABLE WITH STATISTICAL ANALYSIS RESULTS


library(dplyr)
library(gt)
library(gtsummary)
library(WRS2)



# Yuen-Welch t-test
yuen_welch_t_test <- function(data, variable, by, ...) {
  data <- data[c(variable, by)] %>% dplyr::filter(complete.cases(.))
  yuen(data[[variable]] ~ factor(data[[by]]), tr=0.1)$p.value
}

# Permutation Yuen-Welch t-test
permutation_yuen_welch_t_test <- function(data, variable, by, ...) {
  data <- data[c(variable, by)] %>% dplyr::filter(complete.cases(.))
  original_test <- yuen(data[[variable]] ~ factor(data[[by]]), tr=0.1)
  original_statistic <- original_test$test
  n_permutations <- 1000
  permuted_statistics <- numeric(n_permutations)
  
  
  set.seed(123)  # For reproducibility
  for (i in 1:n_permutations) {
    # Permute group labels
    shuffled_group <- sample(data[[by]])
    
    # Calculate the test statistic for the permuted data
    permuted_test <- yuen(data[[variable]] ~ shuffled_group, data = data, tr=0.1)
    permuted_statistics[i] <- permuted_test$test
  }
  
  # Calculate the p-value
  p.value <- mean(abs(permuted_statistics) >= abs(original_statistic))
  
}


# Brown-Mood test of medians
library(RVAideMemoire)
brown_mood_test <- function(data, variable, by, ...) {
  data <- data[c(variable, by)] %>% dplyr::filter(complete.cases(.))
  mood.medtest(data[[variable]] ~ factor(data[[by]]))$p.value
}

# Mean ratios
library(effectsize)
mean_ratios <- function(data, variable, by, ...) {
  data <- data[c(variable, by)] %>% dplyr::filter(complete.cases(.))
  effectsize::means_ratio(data[[variable]] ~ factor(data[[by]]))$Means_ratio
}





table <- tbl_summary(
  data = data%>%
   dplyr::select(-`Gene Name`,-`hsa-miR-208a-3p`, -`hsa-miR-184`,-`hsa-miR-372-3p`,-`hsa-miR-373-3p`, -`hsa-miR-31-5p`, -`hsa-miR-124-3p`, -`hsa-miR-196a-5p`, -`hsa-miR-200b-3p`, -`hsa-miR-203a-3p`, -`hsa-miR-499a-5p`, -`hsa-miR-211-5p`, -`hsa-miR-9-5p`, -`hsa-miR-206`),
     by = Group_level,  # Grouping variable
  type = all_continuous() ~ "continuous2",  # Assuming all biomarkers are continuous
  statistic = all_continuous() ~ c("{N_nonmiss}","{N_miss}","{mean}", "{median} ({p25}-{p75})", "{min} - {max}"), missing="no")%>%
  add_stat_label(label = all_continuous() ~ c("N","N missing", "Mean", "Median (Q1-Q3)", "Min - Max"))%>%
  modify_header(label ~ "**Variable**") %>%
  add_stat(fns = everything() ~ mean_ratios)%>%
  modify_header(add_stat_1 ~ "**Ratio of means**")%>%
  add_p(test= list(all_continuous() ~ permutation_yuen_welch_t_test),
        pvalue_fun = ~ style_pvalue(.x, digits = 3))%>%
  modify_header(p.value ~ "**P-value**")%>%
  modify_footnote(p.value ~ "Permutation Yuen-Welch t test")%>%
  add_q(method = "BH")%>%
  modify_header(q.value ~ "**Multiplicity adjusted p-value**")%>%
  modify_caption("**Table. MicroRNA comparison between study groups**") %>%
  bold_labels()




library(flextable)
table_flextable <- as_flex_table(table)
library(officer)
doc <- read_docx() %>% 
  body_add_flextable(table_flextable)

file_path <- "C:/Users/adawys/.../Report.docx"  
print(doc, target = file_path)









sessionInfo()
# other attached packages:
#  [1] WRS2_1.1-5      gtsummary_1.7.2 gt_0.10.1       rlang_1.1.1     gridExtra_2.3   ggplot2_3.5.0   writexl_1.4.2   dplyr_1.1.2     readxl_1.4.2   
# [10] MatchIt_4.5.5  









