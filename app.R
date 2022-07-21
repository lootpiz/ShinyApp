library(ggplot2)
library(shiny)

tcga_data <- readRDS(file.path("./data", "TCGA_Expression.RDS"))
tcga_gtex <- readRDS(file.path("./data", "GTEx_Expression.RDS"))
tcga_nmfr <- readRDS(file.path("./data", "TCGA_NMF.RDS"))
tcga_rank <- 5

cptac_data <- log2(readRDS(file.path("./data", "CPTAC_Expression.RDS")))
cptac_gtex <- cptac_data[,c(100:108)]
cptac_data <- cptac_data[,c(1:99)]
cptac_nmfr <- readRDS(file.path("./data", "CPTAC_NMF.RDS"))
cptac_rank <- 4

tfri_data <- readRDS(file.path("./data", "TFRI_Expression.RDS"))
tfri_gtex <- NULL
tfri_nmfr <- readRDS(file.path("./data", "TFRI_NMF.RDS"))
tfri_rank <- 3

colours <- c("gold", "orange", "brown1", "brown4", "indianred3")

ui <- fluidPage(
        titlePanel("GBM Bulk RNA-seq"),
        sidebarLayout(
                sidebarPanel(
                        helpText("Query the level of gene expression in either TCGA, CPTAC, or TFRI."),
                        radioButtons(
                                inputId = "dataset", 
                                label = h3("Dataset"),
                                choices = list(
                                        "TCGA" = 1, 
                                        "CPTAC" = 2,
                                        "TFRI" = 3
                                ),
                                selected = 1
                        ),
                        textInput(
                                inputId = "gene",
                                label = h3("Gene symbol"),
                                value = "TNC"
                        ),
                        downloadLink(
                                outputId = "export",
                                label = "Export Figure(s)"
                        )
                ),
                mainPanel(
                        h2(textOutput("selected_dataset")),
                        h5(textOutput("help_dataset")),
                        tabsetPanel(
                                type = "tabs",
                                tabPanel("NMF subgroups", h5("Subgroups were defined via NMF (Non-negatvie Matrix Factorization)"), plotOutput("violin")),
                                tabPanel("Tumour vs. Normal", plotOutput("violin2"))
                        )                        
                )
        )
)

server <- function(input, output, session) {
        generateViolinNMF <- function() {
                if (input$dataset == 1) {
                        data <- tcga_data
                        gtex <- tcga_gtex
                        nmfr <- tcga_nmfr
                        rank <- tcga_rank
                } else if (input$dataset == 2) {
                        data <- cptac_data
                        gtex <- cptac_gtex
                        nmfr <- cptac_nmfr
                        rank <- cptac_rank
                } else if (input$dataset == 3) {
                        data <- tfri_data
                        gtex <- tfri_gtex
                        nmfr <- tfri_nmfr
                        rank <- tfri_rank
                }

                gene <- input$gene

                if (length(which(rownames(data) == gene)) > 0) {
                        if (!is.null(gtex)) {
                                results <- data.frame(
                                        Group = "GTEx",
                                        Expression = as.numeric(gtex[which(rownames(gtex) == gene), ])
                                )
                        } else {
                                results <- c()
                        }

                        subgroups <- nmfr

                        for (subGroupLv in c(1:rank)) {
                                dataSubIdx <- which(colnames(data) %in% names(subgroups)[which(subgroups == subGroupLv)])

                                df <- data.frame(
                                        Group = paste0("Group_", subGroupLv),
                                        Expression = as.numeric(unlist(data[which(rownames(data) == gene), dataSubIdx]))
                                )

                                results <- rbind(results, df)
                        }

                        if (!is.null(gtex)) {
                                lev <- c(paste0("Group_", c(1:rank)), "GTEx")
                                col <- c(colours[1:rank], "grey80")
                        } else {
                                lev <- c(paste0("Group_", c(1:rank)))
                                col <- c(colours[1:rank])
                        }

                        results$Group <- factor(results$Group, levels = lev)

                        dodge <- position_dodge(width = 0.5)
                        ggplot(data = results, aes(x = Group, y = Expression, fill = Group)) +
                                geom_violin(position = dodge, size = 0) +
                                geom_boxplot(width = 0.1, position = dodge, fill = "white") +
                                scale_fill_manual(values = col) +
                                labs(
                                        title = gene,
                                        x = "", y = "Expression"
                                ) +
                                theme_bw() +
                                theme(
                                        axis.line = element_line(colour = "black"),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                                        legend.position = "none",
                                        text = element_text(size = 20)
                                )
                }
        }
        generateViolinTvsB <- function() {
                if (input$dataset == 1) {
                        data <- tcga_data
                        gtex <- tcga_gtex
                        nmfr <- tcga_nmfr
                        rank <- tcga_rank
                } else if (input$dataset == 2) {
                        data <- cptac_data
                        gtex <- cptac_gtex
                        nmfr <- cptac_nmfr
                        rank <- cptac_rank
                } else if (input$dataset == 3) {
                        data <- tfri_data
                        gtex <- tfri_gtex
                        nmfr <- tfri_nmfr
                        rank <- tfri_rank
                }

                gene <- input$gene

                if (length(which(rownames(data) == gene)) > 0) {
                        resultsGbm <- data.frame(
                                Group = "GBM",
                                Expression = as.numeric(data[which(rownames(data) == gene), ])
                        )

                        if (!is.null(gtex)) {        
                                resultsGtex <- data.frame(
                                        Group = "Brain",
                                        Expression = as.numeric(gtex[which(rownames(gtex) == gene), ])
                                )
                                results <- rbind(resultsGbm, resultsGtex)
                                lev <- c("GBM", "Brain")
                                col <- c("chartreuse2", "grey80")
                        } else {
                                results <- resultsGbm
                                lev <- c("GBM")
                                col <- c("chartreuse2")
                        }

                        results$Group <- factor(results$Group, levels = lev)

                        dodge <- position_dodge(width = 0.5)
                        ggplot(data = results, aes(x = Group, y = Expression, fill = Group)) +
                                geom_violin(position = dodge, size = 0) +
                                geom_boxplot(width = 0.1, position = dodge, fill = "white") +
                                scale_fill_manual(values = col) +
                                labs(
                                        title = gene,
                                        x = "", y = "Expression, Log2"
                                ) +
                                theme_bw() +
                                theme(
                                        axis.line = element_line(colour = "black"),
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.border = element_blank(),
                                        panel.background = element_blank(),
                                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                                        legend.position = "none",
                                        text = element_text(size = 20)
                                )
                }
        }
        output$violin <- renderPlot({ generateViolinNMF() }, height = 600, width = 800)

        output$violin2 <- renderPlot({ generateViolinTvsB() }, height = 600, width = 800)

        output$export <- downloadHandler(
                filename = function() {
                        if (input$dataset == 1) {
                                dName <- "TCGA"
                        } else if (input$dataset == 2) {
                                dName <- "CPTAC"
                        } else if (input$dataset == 3) {
                                dName <- "TFRI"
                        } else {
                                dName <- "Wrong"
                        }
                        paste0(dName, "_", input$gene, "_", Sys.Date(), ".pdf")
                },
                content = function(file) {
                        pdf(file = file)
                        print(generateViolinNMF())
                        print(generateViolinTvsB())
                        dev.off()
                }
        )

        output$selected_dataset <- renderText({
                if (input$dataset == 1) {
                        "TCGA-GBM, N=164 vs. GTEx-Cortex, N=105"
                } else if (input$dataset == 2) {
                        "CPTAC-GBM, N=99 vs. GTEx-Brain, N=10"
                } else if (input$dataset == 3) {
                        "TFRI-Tumor, N=44"
                } else {
                        "Something went wrong..."
                }
        })

        output$help_dataset <- renderText({
                if (input$dataset == 1) {
                        "A combined cohort of TCGA and GTEx samples downloaded from Xenabrowser. Gene expression profile is RSEM expected_count output normalized using DESeq2."
                } else if (input$dataset == 2) {
                        "Prospectively collected tumor tissue samples from 99 GBM patients and 10 normal brain samples from CPTAC and GTEx, respectively."
                } else if (input$dataset == 3) {
                        "No (normal) brain profile available. Shen, Grisdale, Islam et al. PNAS (2019). DOI:10.1073/pnas.1813495116"
                } else {
                        "Something went wrong..."
                }
        })
}

shinyApp(ui = ui, server = server)