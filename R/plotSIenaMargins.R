# ## Plotting functions will probably not stay here this way, but for now usable

# ##@plotDiff
# plotDiff <- function(df,
#                     y,
#                     x,
#                     group = NULL,
#                     facet_density = TRUE,
#                     diffvarnames = NULL,
#                     uncertainty = FALSE){
#   if(is.null(diffvarnames)) diffvarnames <- ""
#   if(is.null(group)){
#     p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
#       geom_point() +
#       geom_line() +
#       ylab(paste0("E[",y, "(",paste(diffvarnames, sep = " "),")]"))
#   }else{
#     df[[group]] <- as.factor(df[[group]])
#     p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], 
#                         group = .data[[group]], color = .data[[group]], fill = .data[[group]])) +
#       geom_point() +
#       geom_line() +
#       ylab(paste0("E[",y, "(",paste(diffvarnames, collapse = ","),")]"))
#   }
#   if(uncertainty){
#     # generalize later
#     p <- p +
#       geom_ribbon(aes(ymin= q_025, ymax=q_975), alpha=0.3, linetype="dotted")
#   }
  
#   if(facet_density){
#     p <- p +
#       facet_wrap(~ density, labeller=label_both)
#   }
#   p
# }
# ##@create_heatmap
# create_heatmap <- function(inputmatrix, palette_name = "Blues", scale_title = "probability", summary_func = "mean") {
  
#   prob_df <- as_tibble(inputmatrix, rownames = "ego")
#   prob_df <- gather(prob_df, key = "choice", value = "value", -ego, factor_key = TRUE) %>%
#     mutate(ego = as.numeric(ego), choice = as.numeric(choice))
  
#   prob_df <- bind_rows(prob_df,
#                        summarizeValue(prob_df, choice, ego, summary_func),
#                        summarizeValue(prob_df, ego, choice, summary_func)
#   ) %>%
#     mutate(ego = as.factor(ego), choice = as.factor(choice))
#   prob_heatmap <- ggplot(prob_df, aes(x = choice, y = ego)) +
#     geom_raster(data = filter(prob_df, choice != 0, ego != 0),
#                 aes(fill = value)) +
#     coord_equal() +
#     scale_y_discrete(limits = rev) +
#     scale_fill_distiller(scale_title,
#                          palette = palette_name, direction = 1, limits = c(0, 1)) +
#     theme_void()
#   prob_heatmap +
#     geom_text(data = filter(prob_df, choice == 0 || ego == 0), 
#               aes(label = str_replace((round(value, digits = 3)), "^0|^-0", ""), 
#                   angle = ifelse(ego == 0, -90, 0)), 
#               size = 2)
# }

# ##@create_signed_heatmap
# create_signed_heatmap <- function(inputmatrix,
#                                   palette_name = "RdBu", scale_title = "Tie Change Statistic", summary_func = "mean") {
  
#   prob_df <- as_tibble(inputmatrix, rownames = "ego")
#   prob_df <- gather(prob_df, key = "choice", value = "value", -ego, factor_key = TRUE) %>%
#     mutate(ego = as.numeric(ego), choice = as.numeric(choice))
#   prob_df <- bind_rows(prob_df,
#                        summarizeValue(prob_df, choice, ego, summary_func),
#                        summarizeValue(prob_df, ego, choice, summary_func)
#   ) %>%
#     mutate(ego = as.factor(ego), choice = as.factor(choice))
#   prob_heatmap <- ggplot(prob_df, aes(x = choice, y = ego)) +
#     geom_raster(data = filter(prob_df, choice != 0, ego != 0),
#                 aes(fill=value)) +
#     coord_equal() +
#     scale_y_discrete(limits=rev) +
#     scale_fill_distiller(scale_title,
#                          palette  = palette_name, direction = 1, limits = c(-1,1)) +
#     theme_void()
#   prob_heatmap +
#     geom_text(data = filter(prob_df, choice == 0 || ego == 0),
#               aes(label = str_replace((round(value, digits = 3)), "^0|^-0", ""),
#                   color = (value>=0),
#                   angle = ifelse(ego == 0, -90, 0)), 
#               size = 2) +
#     scale_colour_manual(summary_func,
#                         values = RColorBrewer ::brewer.pal(11,palette_name)[c(3, 9)],
#                         guide = FALSE)
# }

# ##@save_plot
# save_plot <- function(rel_path, ggp, width = 60, height = 36, dpi = 300, scaling = 3){
#   width_px <- width * 0.3937 * dpi
#   height_px <- height * 0.3937 * dpi
#   res <- dpi * scaling
#   ggsave(rel_path, ggp, width = width_px, height = height_px, units = "px", dpi = res)
# }
  
