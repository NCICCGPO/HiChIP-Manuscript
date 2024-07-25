#get heatmap color gradient from any color ramp to white, with specified dark portion for emphasis
getColorRampGradient <- function(dark_color = NA, 
                                 mid_color = NA, 
                                 viridis_palette = 'magma', 
                                 light_color = 'white', 
                                 paletteLength = 50, 
                                 color_dark_portion = 2/3){ 
    #top 2/3 of values will have darker colors on the heatmap unless color_dark_portion is specified by user
    paletteLength_dark_portion <- round( paletteLength*color_dark_portion )
    paletteLength_light_portion <- paletteLength - paletteLength_dark_portion
    colors_were_provided <- c(!is.na(dark_color), !is.na(mid_color))
    if ( all(colors_were_provided) ){ 
        dark_to_mid <- colorRampPalette( c(dark_color,mid_color) )(paletteLength_dark_portion)
    } else {
        require(viridis)
        dark_to_mid_dictionary <- list('magma' = magma(paletteLength_dark_portion),
                                       'cividis' = cividis(paletteLength_dark_portion),
                                       'viridis' = viridis(paletteLength_dark_portion),
                                       'inferno' = inferno(paletteLength_dark_portion),
                                       'plasma' = plasma(paletteLength_dark_portion))
        dark_to_mid <- dark_to_mid_dictionary[[viridis_palette]]
    }
    mid_to_light <- colorRampPalette(c(dark_to_mid[length(dark_to_mid)], light_color))(paletteLength_light_portion)
    mycolors <- rev(c(dark_to_mid, mid_to_light))
    return(mycolors)
}

#generate color function using colorRamp2
generate_colorRamp2_colfun <- function(vector, 
                                       dark_color_portion = 1, 
                                       min_value = min(vector),
                                       mid_value = max(vector) - (max(vector) - min(vector)) * dark_color_portion,
                                       max_value = max(vector),
                                       min_mid_max_colors = c("white", "white", "firebrick")){
    min_mid_max <- c(min_value, mid_value, max_value)
    col_fun <- circlize::colorRamp2(min_mid_max, min_mid_max_colors)
    return(col_fun)
}

#custom palettes
splash_of_salmon <- c('#e05c4d', '#bad532', '#7c953b', '#154d93', '#477cbf')
yennefer_violet <- c('#e66101', '#fdb863', '#f7f7f7', '#b2abd2', '#5e3c99')
crookback_bog <- c('#486773', '#7297a6', '#d9c6b0', '#d99066', '#bf6e50')
firebrick_divergent <- c("navy", "white", "firebrick")
blue_orange_divergent <- c("#3794bf", "#FFFFFF", "#df8640")

my_custom_palettes <- list(splash_of_salmon = splash_of_salmon, 
                           yennefer_violet = yennefer_violet, 
                           crookback_bog = crookback_bog,
                           firebrick_divergent = firebrick_divergent,
                           blue_orange_divergent = blue_orange_divergent)