library(dplyr)
library(tidyr)
library(ggcyto)
library(ggplot2)
library(stringr)
library(magrittr)
library(flowCore)
library(openCyto)



## CONFIG B  Before 10/20, Config A After 10/20
flurophore_marker_mapping <-
    c(
        "PE-Cy7"          = "CD13",
        "PeCy7"           = "CD13",
        "FITC"            = "CD15",
        "BB515"           = "CD15",
        "BV605"           = "CD19",
        "BV421"           = "CD19",
        "PE"              = "CD33",
        "APC"             = "CD34",
        "APC-ALEXA 750"   = "CD38",
        "APC-H7"          = "CD38",
        "APC-Cy7"         = "CD38",
        "V500"            = "CD45",
        "BUV805"          = "CD45",
        "Alexa 700"       = "CD71",
        "Alexa Fluor 700" = "CD71",
        "PerCP-Cy5"       = "CD117",
        "PE-Cy5"          = "CD117",
        "BB700"           = "CD117",
        "Pacific Blue"    = "HLADR",
        "PACBLU"          = "HLADR",
        "BUV395"          = "HLADR",
        "BUV396"          = "HLADR"
    )



canto_lot_e_targets <-
    c(
        "FITC"            = 37045,
        "PE"              = 72865,
        "PerCP-Cy5"       = 48787,
        "PeCy7"           = 11282,
        "APC"             = 128165,
        "ALEXA FLUOR 700" = 64706,
        "APC-ALEXA 750"   = 22616,
        "PACIFIC BLUE"    = 89652,
        "V500"            = 139416,
        "BV605"           = 46791
    )

canto_lot_a_targets <-
    c(
        "FITC"            = 31876,
        "PE"              = 68495,
        "PerCP-Cy5"       = 42983,
        "PeCy7"           = 9880,
        "APC"             = 86745,
        "ALEXA FLUOR 700" = 40107,
        "APC-ALEXA 750"   = 14032,
        "PACIFIC BLUE"    = 148216,
        "V500"            = 158977,
        "BV605"           = 63303
    )


fortessa_before_9_20_lot_a_targets <-
    c(
        "FITC"            = 65460,
        "PE"              = 103616,
        "PerCP-Cy5"       = 43464,
        "PeCy7"           = 14242,
        "APC"             = 103795,
        "ALEXA FLUOR 700" = 29909,
        "APC-ALEXA 750"   = 15554,
        "PACBLU"          = 101369,
        "V500"            = 160368,
        "BV605"           = 95919
    )

# HLA-DR PACBLU PREVIOUS NOW HLA-DR BUV395
# CD45 V500     PREVIOUS NOW CD45 BUV805
# CD5 PercP-Cy5.5 PREVIOUSLY NOW CD5 BB700
# 6TH peak instead of 7th peak
# NOTE: Since we will use target peak values to locate the
# peak, we do not really need to keep track of 6th/7th "peak" use
# information.
fortessa_after_9_20_lot_a_targets <-
    c(
        "FITC"                 = 9632,
        "PE"                   = 24313,
        "PE-Texas Red-A"    = 22262,
        "PE-Cy5-A"          = 11952,
        "PE-Cy7-A"          = 2887,
        "BV421-A"           = 23145,
        "BV510-A"           = 13105,
        "BV605-A"           = 8471,
        "BV650-A"           = 3371,
        "BV711-A"           = 1519,
        "BV785-A"           = 373,
        "APC-A"             = 26288,
        "Alexa Fluor 700-A" = 15263,
        "APC-H7-A"          = 3932,
        "BUV395-A"          = 8366,
        "BUV563-A"          = 7846,
        "BUV737-A"          = 456
    )


symphony_before_10_20_lot_a_targets <-
    c(
        "FITC"            = 65460,
        "PE"              = 103616,
        "PerCP-Cy5"     = 43464,
        "PeCy7"           = 14242,
        "APC"             = 103795,
        "ALEXA FLUOR 700" = 29909,
        "APC-ALEXA 750"    = 15554,
        "PACBLU"          = 101369,
        "V500"            = 160368,
        "BV605"           = 95919
    )

symphony_after_10_20_lot_a_targets <-
    c(
        "BB515-A "     = 4157,
        "BB630-A "     = 18386,
        "BB660-A "     = 9500,
        "BB700-A "     = 3764,
        "BB790-A "     = 3140,
        "APC-A "       = 28350,
        "Alexa 700-A " = 18930,
        "APC-Cy7-A "   = 4119,
        "BUV396-A "    = 8763,
        "BUV496-A "    = 21255,
        "BUV563-A "    = 8579,
        "BUV615-A "    = 7596,
        "BUV661-A "    = 4150,
        "BUV737-A "    = 559,
        "BUV805-A "    = 161,
        "BV421-A "     = 15788,
        "BV480-A "     = 20770,
        "BV570-A "     = 10698,
        "BV605-A "     = 5036,
        "BV650-A "     = 2952,
        "BV711-A "     = 1646,
        "BV750-A "     = 883,
        "BV786-A "     = 372,
        "PE-A "        = 26139,
        "PE-CF594-A "  = 29372,
        "PE-Cy5-A "    = 16725,
        "PE-Cy5-5-A "  = 7064,
        "PE-Cy7-A "    = 4384
    )





map_fluoros_to_proteins <- function(df) {

    flurophore_marker_mapping <-
        flurophore_marker_mapping[
            order(str_length(names(flurophore_marker_mapping)))
        ]

    protein_factor_levels <-
        c("CD13", "CD15", "CD19", "CD33", "CD34", "CD38", "CD45",
        "CD71", "CD117", "HLADR")

    for (fluorophore in names(flurophore_marker_mapping)) {

        reg_pattern <- regex(fluorophore, ignore_case = TRUE)
        detected <- str_detect(df$name, reg_pattern)

        if (!any(detected)) next

        channels <- df$name[detected]
        lengths <- str_length(channels)
        channel <- channels[lengths == min(lengths)]


        df[df$name == channel, 'protein'] <-
            flurophore_marker_mapping[fluorophore]

    }

    df <- df %>%
        tidyr::drop_na() %>%
        mutate(protein = factor(protein, protein_factor_levels)) %>%
        arrange(protein)

    return(df)
}




read_fcs_data <- function(fcs_files) {

    n_files <- length(fcs_files)
    fcs_data <- vector('list', n_files)

    for (i in seq_along(fcs_files)) {
        cat('\nReading FCS file: ', i, ' of: ', n_files, '\n')

        fcs_data[[i]] <-
            suppressMessages(
                read.FCS(fcs_files[[i]],
                    transformation = FALSE,
                    truncate_max_range = FALSE))
    }

    return(fcs_data)
}




gate_singlets <-
    function(fcs,
              plot_gate = FALSE) {

    FF <- flowFrame(exprs(fcs))

    # Gate for singlets and viable
    chnl <- c("FSC-A", "FSC-H")
    singlets_gate <-
        suppressMessages(openCyto:::.singletGate(FF,
                                channels = chnl,
                                wider_gate = TRUE))

    if (isTRUE(plot_gate)) {
        p <- autoplot(FF,
                    x = chnl[1],
                    y = chnl[2],
                    bins = 100)

        p1 <- p + geom_gate(singlets_gate)
        p1 <- as.ggplot(p1)
        print(p1)
    }

    # Selecting Singlets and continuing
    singlets_filter <-
        suppressMessages(
            flowCore::filter(FF, singlets_gate))

    singlets <-
        suppressMessages(
            flowCore::split(FF, singlets_filter)[[1]]
            )

    prop <- nrow(exprs(singlets)) / nrow(exprs(fcs))
    cat('Proportion singlets selected: ', prop, '\n')

    return(singlets)
}







gate_good_beads <-
    function(fcs,
             plot_gate = FALSE) {

    flowclust_gate <-
        suppressMessages(
            openCyto:::.flowClust.2d(
                fcs,
                channels = c('FSC-A', 'SSC-A'),
                K=1,
                target=c(70e3,5e4), quantile=0.95))

    filter_from_gate <-
        suppressMessages(
            flowCore::filter(fcs, flowclust_gate))

    good_beads  <-
        suppressMessages(
            flowCore::split(fcs, filter_from_gate) %>%
            .[[1]])

    if (isTRUE(plot_gate)) {

        new_exprs <- exprs(fcs)
        new_exprs[new_exprs <= 0] = 0.1
        plot_fcs <- flowFrame(exprs = new_exprs)

        p <-
            autoplot(
                plot_fcs,
                x = 'FSC-A',
                y = 'SSC-A',
                bins = 100) +
                scale_x_log10() +
                scale_y_log10()

        p <-
            p + geom_gate(flowclust_gate)

        print(p)

    }

    return(good_beads)
}





gate_target_population <-
    function(fcs,
             peak_targets = c('PE-A'   = 70e3,
                              'FITC-A' = 35e3),
             plot_gate = FALSE) {


   if (!all(names(peak_targets) == c('PE-A', 'FITC-A'))) {
       stop(paste0('\nERROR: peak_targets must have names PE-A, FITC-A, not: ', names(peak_targets)))
   }

    p7_gate <-
        suppressMessages(
            openCyto:::.flowClust.2d(
                fcs,
                channels = c('PE-A', 'FITC-A'),
                K = 8,
                target = peak_targets,
                quantile = 0.99
            )
        )

    if (isTRUE(plot_gate)) {

        new_exprs <- exprs(fcs)
        new_exprs[new_exprs <= 0] = 0.1
        plot_fcs <- flowFrame(exprs = new_exprs)

        p <-
            autoplot(
                plot_fcs,
                x = 'PE-A',
                y = 'FITC-A',
                bins = 100) +
            geom_gate(p7_gate) +
            geom_point(aes(x = peak_targets[1],
                           y = peak_targets[2]),
                           shape = 23, size = 2) +
            xlim(c(peak_targets[1] * .1, peak_targets[1] * 2)) +
            ylim(c(peak_targets[2] * .1, peak_targets[2] * 2)) +
            scale_x_continuous(trans = 'log10', labels = scales::comma) +
            scale_y_continuous(trans = 'log10', labels = scales::comma)


        print(p)
    }

    # Selecting gate and continuing
    filter_from_gate <-
        suppressMessages(
            flowCore::filter(fcs, p7_gate)
        )

    p7_cells <-
        suppressMessages(
            flowCore::split(fcs, filter_from_gate) %>%
            .[[1]])

    return(p7_cells)
}









find_fcs_files_from_qc_dirs <- function(qc_dirs) {

    n_dirs <- length(qc_dirs)
    fcs_files <- vector('list', n_dirs)

    for (i in seq_along(qc_dirs)) {

        cat('Getting FCS files for dir: ', i, ', of: ', n_dirs, '\n')

        sphero_qc_subdir <-
            qc_dirs[[i]] %>%
            list.dirs(full.names = TRUE) %>%
            .[stringr::str_detect(., regex('daily qc', ignore_case = TRUE))]

        if (length(sphero_qc_subdir) == 1) {
            fcs_file_to_use <- list.files(sphero_qc_subdir,
                                         pattern = '\\.fcs',
                                         full.names = TRUE)

            if (length(fcs_file_to_use) > 1) {
                file_mtimes <- file.info(fcs_file_to_use)$mtime
                last_file <- which.max(file.info(fcs_file_to_use)$mtime %>% as.integer())
                fcs_file_to_use <- fcs_file_to_use[last_file]
            }

             fcs_files[[i]] <- fcs_file_to_use

        } else {
            message(paste0('Issue found with directory: ', qc_dirs[[i]]))
        }
    }

    return(fcs_files)
}




retrieve_peak_means <- function(fcs_data,
                                peak_targets = c('PE-A' = 72865, 'FITC-A' = 37045),
                                plot_gate = FALSE) {

    singlets   <- lapply(fcs_data, function(fcs) gate_singlets(fcs))
    good_beads <- lapply(singlets, function(fcs) gate_good_beads(fcs))
    target_cells <- lapply(good_beads,
        function(fcs) {
            gate_target_population(
                fcs,
                plot_gate = plot_gate,
                peak_targets = peak_targets)
        })

    peak_means <- lapply(target_cells,
        function(fcs) {
             as.data.frame(t(colMeans(exprs(fcs))))
        })

    peak_means <- do.call(rbind, peak_means)

    return(peak_means)
}





# Example Workflow

if (FALSE) {

    qc_dirs <- c('YOUR QC DIRS HERE')
    fcs_files <- find_fcs_files_from_qc_dirs(qc_dirs)
    fcs_data <- read_fcs_data(fcs_files)
    singlets   <- lapply(fcs_data, function(fcs) gate_singlets(fcs))
    good_beads <- lapply(singlets, function(fcs) gate_good_beads(fcs))
    target_cells <- lapply(good_beads, function(fcs) gate_target_population(fcs, plot_gate = TRUE, peak_targets = c('PE-A' = 24313, 'FITC-A' = 9632)))

    peak_means <- lapply(target_cells, function(fcs) as.data.frame(t(colMeans(exprs(fcs)))))
    peak_means <- do.call(rbind, peak_means)
    peak_means


    canto_run_data <- list(canto_fcs_1, canto_fcs_2, canto_fcs_3)
    canto_peak_means <- retrieve_peak_means(canto_run_data, peak_targets = c('PE-A' = 72865, 'FITC-A' = 37045), TRUE)
    canto_peak_means %>% select(`FITC-A`, `PE-A`)

    fortessa_run_data <- list(fortesaa_fcs_1, fortesaa_fcs_2, fortesaa_fcs_3)
    fortessa_peak_means <- retrieve_peak_means(fortessa_run_data, peak_targets = c('PE-A' = 24313, 'FITC-A' = 9632), TRUE)
    fortessa_peak_means %>% select(`FITC-A`, `PE-A`)

}
