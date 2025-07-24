#' @export
lzq_venn.plot <- function (data.list, set.names = NULL, fill.colors = c("#BC3C29FF", 
    "#0072B5FF", "#E18727FF", "#20854EFF"), fill.alpha = 0.5, 
    border.color = "black", border.alpha = 1, border.width = 1, 
    border.linetype = 1, set.name.colors = "black", set.name.size = 6, 
    text.color = "black", text.size = 5.5) 
{
    if (length(data.list) <= 1) {
        stop("The number of elements in data.list must be greater than 1!")
    }
    if (length(data.list) > 5) {
        stop("This function currently only supports venn plot of classes 2-4!")
    }
    if (is.null(set.names)) {
        if (is.null(data.list)) {
            stop("Elements in data.list must have names, or set.names must be prepared!")
        }
        else {
            set.names <- names(data.list)
        }
    }
    data.length <- length(set.names)
    data <- purrr::map2_df(data.list, names(data.list), ~data.frame(value = .x, 
        Class = .y, Type = TRUE)) %>% tidyr::pivot_wider(names_from = "Class", 
        values_from = "Type") %>% dplyr::mutate(dplyr::across(dplyr::everything(), 
        ~tidyr::replace_na(., FALSE)))
    venn_plot <- geom_venn(set_names = set.names, fill_color = fill.colors, 
        fill_alpha = fill.alpha, stroke_color = border.color, 
        stroke_alpha = border.alpha, stroke_linetype = border.linetype, 
        stroke_size = border.width, set_name_color = set.name.colors, 
        set_name_size = set.name.size, text_color = text.color, 
        text_size = text.size, show_percentage = FALSE, digits = 1, 
        show_outside = "none", )
    if (data.length == 2) {
        colnames(data)[-1] <- LETTERS[1:2]
        p <- ggplot(data, aes(A = A, B = B)) + theme_void() + 
            coord_fixed() + venn_plot
    }
    if (data.length == 3) {
        colnames(data)[-1] <- LETTERS[1:3]
        p <- ggplot(data, aes(A = A, B = B, C = C)) + theme_void() + 
            coord_fixed() + venn_plot
    }
    if (data.length == 4) {
        colnames(data)[-1] <- LETTERS[1:4]
        p <- ggplot(data, aes(A = A, B = B, C = C, D = D)) + 
            theme_void() + coord_fixed() + venn_plot
    }
    print(p)
}
